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
# Function creates a text file with all project details
# This function is based on an original Matlab function development by 
# A. A. Bloom (UoE, now at the Jet Propulsion Laboratory). Translation to R and 
# subsequent modifications by T. L Smallman (t.l.smallman@ed.ac.uk, UoE).
#
#########################################################################################

cardamom_write_project_info <- function (PROJECT) {

    # create file name and output basic information
    file=paste(PROJECT$localpath,"/",PROJECT$name,"_run_DETAILS.txt", sep="")
    write(paste("Name: ",PROJECT$name,sep=""),sep=",",ncolumns = 1,file=file,append="F")
    write(paste("Number of sites: ",PROJECT$nosites,sep=""),sep=",",ncolumns = 1,file=file,append="T")
    write(paste("Number of MCMC chains: ",PROJECT$nochains,sep=""),sep=",",ncolumns = 1,file=file,append="T")
    write(paste("Number of accepted samples: ",PROJECT$nsamples,sep=""),sep=",",ncolumns = 1,file=file,append="T")
    write(paste("Number of accepted subsamples: ",PROJECT$nsubsamples,sep=""),sep=",",ncolumns = 1,file=file,append="T")
    write(paste("Description: ",PROJECT$description,sep=""),sep=",",ncolumns = 1,file=file,append="T")
    write(paste("Date: ",PROJECT$date,sep=""),sep=",",ncolumns = 1,file=file,append="T")

    # load and output parameter information
    pdetails = parameter_details(PROJECT$model$name,PROJECT$parameter_type,PROJECT$ctessel_pft)
    ranges = parameter_ranges(PROJECT$model$name,PROJECT$parameter_type,PROJECT$ctessel_pft)
    Pmin = ranges$allparsmin ; Pmax = ranges$allparsmax
    for (n in seq(1, length(Pmin))) {
         write(paste(pdetails$parameter_names[n],"min and max values ",Pmin[n]," - ",Pmax[n],sep=""),sep=",",ncolumns = 1,file=file,append="T")
    }
}

## Use byte compile
cardamom_write_project_info<-cmpfun(cardamom_write_project_info)
