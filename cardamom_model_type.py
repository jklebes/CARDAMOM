# 1) insert 'type model_working_variables' and 'end type' around the 
#    non-parameter module-level variables collection
# 2) Get rid of all the variables in "public:" block except subroutiens and parameters: 
#    CARBON_MODEL and nos_soil_layers, top_soil_depth, other read-only parameters USEd by likelihood file
# 3) run this script
#    'python cardmom_model_type.py filename.f90'
#
#       - changes references  from x to mv%x
#       - inserts mV into subroutine args list where needed
#       - inserts mV declaration into subroutines where needed
#
#     (doesn't overwrite, writes to _editted file )
#
# 4) Copy altered file to original filename
# 5) Correct function calls that were missed


import re
import sys

# read file
with open(sys.argv[1], 'r') as f:
    lines= f.readlines()
f.close()

# find start and end of 'model_working_variables' type
started = False
for i,line in enumerate(lines):
    if "type model_working_variables" in line:
        var_list_start = i
        started = True
    if started and "end type" in line:
        var_list_end = i
        started = False
        break
var_lines = lines[var_list_start+1:var_list_end]

# drop comments
def drop_comment(s):
    result = ''
    i = 0
    c = s[i]
    while c != '!' and i < len(s):
        result = result + c
        i+=1
        try:
            c = s[i]
        except:
            pass
    return result
var_lines = [drop_comment(s) for s in var_lines]

# split at whitespace and flatten into one list of words
var_lines = [l.split() for l in var_lines]
def flatten(xss):
    return [x for xs in xss for x in xs]

var_lines = flatten(var_lines)
var_lines_tmp = var_lines 
var_lines = []
for word in var_lines_tmp:
    if len(word)>=2 and word[-2:]=='::':
        var_lines.append(word[:-2])
        var_lines.append(word[-2:])
    elif len(word)>=1 and word[-1:]==',':
        var_lines.append(word[:-1])
        var_lines.append(word[-1:])
    else:
        var_lines.append(word)

# keep variable names after '::',  ',' or  ', &'
var_section = False
var_names = []
word_count = 0
for word in var_lines:
    if ',' in word or '&' in word:
        word_count = 0
    else:
        word_count +=1
    if word_count > 1 :
        var_section = False
    if var_section :
        if ',' not in word and '&' not in word and word !="double":
            var_names.append(word)
            word_count  += 1
    if word == '::':
        var_section = True
        word_count = 0

print("Found type variables to change: ", var_names)

# rewrite file:
# after type declaration: 
# transform every instance of those variables from x to mV%x
def alter_vars(lines, line_start, line_end, var_names):
    for v in var_names:
        expr = re.compile(r"(?i)\b("+v+")\\b")
        for i,line in enumerate(lines[line_start:line_end]):
            # need regex for v, case insensitive, bordered by whitespace or punctuation
            # replace with same prefixed with 'mV%'
            line_p = line.partition('!')
            lines[line_start+i] = re.sub(pattern=expr, repl='mV%\\1' , string=line_p[0]) + line_p[1] + line_p[2]
    return lines

def alter_subroutine(line_start, line_end, use=False):
    if any("mV%" in line for line in lines[line_start+1:line_end]):
        if not use:  
            # if there is an existing args list insert mV before  )
            line = lines[line_start]
            if ')' in line.partition('!')[0]:
                index = line.find(')')
                lines[line_start] = line[:index ]+ ", mV" + line[index:]
            # else insert (mV)
            else:
                subroutine_name = line.strip().split()[1]
                index = line.find(subroutine_name)
                lines[line_start] = line[:index+len(subroutine_name)]+ " (mV)" + line[index+len(subroutine_name):]
        print("altered procedure declaration", lines[line_start])
        declared= False
        # check if there is an `implicit none` in the subroutine body
        for i,l in enumerate(lines[line_start+1:line_end]):
            if "implicit none" in l.partition('!')[0]:
            # insert mV declaration after implicit none
                insert_at = line_start+1+i
                declared = True
        if not declared:
            # else insert mV declaration after any 'use' lines
            for i,l in enumerate(lines[line_start+1:line_end]):
                if not ("use" in l.partition('!')[0] or len(l.partition('!')[0].strip())==0):
                    insert_at = line_start+1+i
                    break
        return insert_at
    return None

def find_subroutine_names(lines):
    # find subroutine starts, ends, vars
    started = False
    collecting_args = False
    subroutine_names = []
    subroutine_starts = []
    subroutine_ends = []
    subroutine_vars = []
    subroutine_type = []
    for i,line in enumerate(lines):
        if not started and any(x in line.partition('!')[0] for x in ("function","subroutine")):
            if "function" in line.partition("!")[0]:
                subroutine_type = "function"
            else:
                subroutine_type = "subroutine"
            subroutine_starts.append(i)
            # read subroutine name and variable names from the line(s)
            l = re.split('\W+', line.strip())
            subroutine_names.append(l[l.index(subroutine_type)+1])
            collecting_args = True
            args = []
            if '(' not in line.partition('!')[0]:
                collecting_args = False
                subroutine_vars.append(args)
            else:
                args_section =  line.partition('!')[0].partition('(')[2].strip() # line before comment and after ( 
                args =  re.split('\W+', args_section)
                if len(args_section)>0 and args_section[-1] == '&':
                    collecting_args = True
                elif len(args_section)>0 and args_section[-1] == ')':
                    collecting_args = False # encountered end of subroutine args list
                    subroutine_vars.append(args)
            started = True
        if started and collecting_args: #in multiline subroutine args list
            subroutine_starts[-1] = i
            args_section =  line.partition('!')[0].strip()
            args = args + re.split('\W+', args_section)
            if len(args_section)>0 and args_section[-1] == '&':
                collecting_args = True
            elif len(args_section)>0 and args_section[-1] == ')':
                collecting_args = False # encountered end of subroutine args list
                subroutine_vars.append(args)
        if started and "end "+subroutine_type in line.partition('!')[0]:
            subroutine_ends.append(i)
            started = False

    print("Found", len(subroutine_names)," procedures:", subroutine_names)
    return (subroutine_names , subroutine_starts, subroutine_ends, subroutine_vars)

def find_insertions(lines, subroutine_info):
    """
    Identify subroutines from the list which use variables that should now be mV, 
    change those variables, 
    report points to insert declarations of mV later 
    """
    subroutine_names, subroutine_starts, subroutine_ends, subroutine_vars = subroutine_info
    insert_declarations = []
    procedures_altered = []
    for i, subroutine_name in enumerate(subroutine_names):
        var_names_subroutine = [v for v in var_names if v not in subroutine_vars[i]]
        lines = alter_vars(lines, subroutine_starts[i]+1, subroutine_ends[i], var_names_subroutine)
        declaration_inserted = alter_subroutine(subroutine_starts[i], subroutine_ends[i])
        if declaration_inserted is not None:
            insert_declarations.append(declaration_inserted) # record line numbers to insert at
            procedures_altered.append(subroutine_name) # record names of procedures now using mV
    return(lines, insert_declarations, procedures_altered)

def edit_original_file(lines, procedures_altered, insert_declarations, filename):
    # edit function calls
    opened = 0
    closed = 0
    cont = False
    for subroutine_name in procedures_altered:
        for i,line in enumerate(lines):
            line_p = line.partition('!')[0]
            words = re.split('\W+', line_p.strip())
            if ("call "+subroutine_name in line_p and "call "+ subroutine_name+"(" not in line_p and 
                 "call "+ subroutine_name+"_" not in line_p ):
                index = line_p.find(subroutine_name) + len(subroutine_name)
                lines[i] = line[:index] + "(mV)" + line[index:]
                print("altered procedure call", lines[i])
            elif cont or "call "+ subroutine_name + "(" in line_p:
                for j,c in enumerate(line_p[max(0,line_p.find(subroutine_name)):]):
                    if c=='(':
                        opened +=1 
                    if c == ')':
                        closed += 1
                        if opened == closed:
                            index = max(0,line_p.find(subroutine_name))+j
                            cont = False
                            opened = 0
                            closed = 0
                            lines[i] = line[:index] + ", mV" + line[index:]
                            print("altered procedure call", lines[i])
                    if c == "&":
                        cont = True

    declaration= "\n      type(model_working_variables) :: mV\n"
    insertions = 0
    lines0 = lines
    lines = []
    # insert declarations into list of lines (shifts line numbering)
    for i,line in enumerate(lines0):
        lines.append(lines0[i])
        if i in insert_declarations:
            lines.append(declaration)
            print("inserting mV declaration at" , i+insertions)
            insertions+=2

    with open(filename+"_editted" , 'w') as f:
        f.writelines(lines)
    f.close()

procedures_altered_all = []
for i,filename in enumerate(sys.argv[1:]):
    with open(filename, 'r') as f:
        lines= f.readlines()
    f.close()
    subroutine_info = find_subroutine_names(lines)
    print("Changing references in file", filename)
    lines, insertion_indices, procedures_altered = find_insertions(lines, subroutine_info)
    procedures_altered_all.extend(procedures_altered)
    if i == 0:
        edit_original_file(lines, procedures_altered_all, insertion_indices, filename)
    else:
        print("Editing additional files with relevant use statemtens not implmented yet")
