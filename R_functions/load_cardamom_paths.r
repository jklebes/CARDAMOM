#########################################################################################
# CARbon DAta MOdel fraMework (CARDAMOM) and DALEC terrestrial ecosystem model suite
# CARDAMOM is a Bayesian model-data fusion software framework. CARDAMOM is used to 
# assimilate observations and ecological theory to retrieve parameters for the 
# DALEC suite of intermediate complexity terrestrial ecosystem models. DALEC can be
# used as a fully integrated component of CARDAMOM or independently. 
# Copyright (C) 2024  University of Edinburgh,
#                     Mathew Williams (mat.williams@ed.ac.uk), 
#                     T. Luke Smallman (t.l.smallman@ed.ac.uk)
# UoE = University of Edinburgh

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

# ########## File specific description ##########
# Loads the relevant file paths for CARDAMOM
# 
# This function is based on an original Matlab function development by A. A. Bloom 
# (UoE, now at the Jet Propulsion Laboratory). Translation to R and subsequent 
# modifications by T. L Smallman (t.l.smallman@ed.ac.uk, UoE).
#
#########################################################################################

load_paths<- function() {

    # cardamom common paths file name
    cardamompathfile="./CARDAMOM_LOCAL/cardamompathdef.RData"
    if(file.exists("./CARDAMOM_LOCAL") == FALSE) {system("mkdir ./CARDAMOM_LOCAL/")}

    if (file.exists(cardamompathfile) == FALSE) {
        # ask some information
        outputsdir=readline("Enter the output location for all your CARDAMOM outputs (e.g. /yourlocaldisk/CARDAMOM/CARDAMOM_OUTPUTS/)")
        cluster=readline("Enter the remote cluster address (e.g. eddie3.ecdf.ed.ac.uk)")
        ecdfdir=readline("Enter the CARDAMOM directory on the remote cluster (e.g. /exports/work/geos_gc_ctessel/CARDAMOM/)")
        # force some slashes
        outputsdir=paste(outputsdir,"/",sep="")
        ecdfdir=paste(ecdfdir,"/",sep="")
        # pretty much all this does
        load_paths=list(user=username,cardamom=paste(getwd(),"/",sep="")
                       ,cardamom_library=paste(getwd(),"/LIBRARY/",sep="")
                       ,cardamom_projects=paste(getwd(),"/PROJECTS/",sep="")
                       ,cardamom_outputs=outputsdir
#                       ,cardamom_data=datadir
                       ,cardamom_cluster=cluster
                      ,cardamom_ecdf=ecdfdir)
        save(load_paths,file=cardamompathfile)
    } else {
        # if file already exists then load from file
        load(cardamompathfile)
    }
    # create directories if needed
    if (file.exists(load_paths$cardamom_outputs) == FALSE ){system(paste("mkdir ",load_paths$cardamom_outputs,sep=""))}
    # output
    return(load_paths)

} # end function load_paths
