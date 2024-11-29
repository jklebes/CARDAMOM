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
# Function automates accessing the UoE HPC service Eddie
# and issues commands to Eddie regarding CARDAMOM
# Could be generically applied to other remote servers with similar security
# Author: T. Luke Smallman (02/05/2024)
#
#########################################################################################

#ecdf_execute<- function(command,cluster_address) {
#
#    # create fifo  connection for ssh
#    if (file.exists("sshOut") == TRUE ){ system('rm sshOut')}
#    system('mkfifo sshOut')
#    #Then, connect to the pipe and the fifo from within R:
#    # The &> redirects both stdout and stderr to the fifo
#    #sshIn <- pipe( paste("ssh -A -tt ",username,"@",cluster_address," &> sshOut",sep=""), open = 'w+')
#    sshIn <- pipe( paste("ssh -A -T ",username,"@",cluster_address,sep=""), open = 'w+')
#    sshOut <- fifo( 'sshOut', 'r+' )
#    # write commands to list
#    for (i in seq(1, length(command))) {writeLines( command[i], sshIn )}
#    # 'flush' pass arguements and then read remote systems return
#    flush( sshIn )
#    readLines( sshOut )
#    # clean up
#    close(sshIn)
#    system('rm sshOut')
#
#} # end function ecdf_execute

ecdf_execute<- function(command,cluster_address) {

    # Clear any existing sshOut file
    if (file.exists("sshIn.sh") == TRUE ){ system('rm sshIn.sh')} 

    # Create the new bash script that contains the commands we want to run
    col_sep = "" ; nos_cols = 2
    write(c("#!/bin/sh"), file = "sshIn.sh", ncolumns = nos_cols, sep=col_sep, append = FALSE)
    # write commands to file
    for (i in seq(1, length(command))) {
         write(command[i], file = "sshIn.sh", ncolumns = nos_cols, sep=col_sep, append = TRUE)
    }
    # Make sure the user is there befor attempting to connect to avoid timing out issues
    here = "n"
    while (here != "y") {
       here = readline("...is there anybody out there...? (y/n)")
    }
    # Issue commands to the remote server
    system(paste("cat sshIn.sh | ssh -i ",sshpass_key_server," ",username,"@",cluster_address,sep=""))
    # Tidy up the command file
    system('rm sshIn.sh')

} # end function ecdf_execute

