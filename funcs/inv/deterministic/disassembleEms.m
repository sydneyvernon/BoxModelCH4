%%% =======================================================================
%%% = disassembleEms.m
%%% = Alex Turner
%%% = 05/03/2017
%%% =----------------------------------------------------------------------
%%% = NOTES
%%% =  ( 1): Break the emissions into a structure.
%%% =  ( 2): "assembleEms" will do the opposite.
%%% =----------------------------------------------------------------------
%%% = INPUTS
%%% =  ( 1): in -- Matrix with the emissions.
%%% =----------------------------------------------------------------------
%%% = OUTPUTS
%%% =  ( 1): out -- Structure with the emissions.
%%% =======================================================================

function [ out ] = disassembleEms( in )
out.nh_ch4   = in(:,1);
out.sh_ch4   = in(:,2);
out.nh_13ch4 = in(:,3);
out.sh_13ch4 = in(:,4);
out.nh_14ch4 = in(:,5);
out.sh_14ch4 = in(:,6);
out.nh_ch3d  = in(:,7);
out.sh_ch3d  = in(:,8);
out.nh_oh    = in(:,9);
out.sh_oh    = in(:,10);
out.nh_co    = in(:,11);
out.sh_co    = in(:,12);
out.tau_TS   = in(:,13);
out.kX_NH    = in(:,14);
out.kX_SH    = in(:,15);
out.cl_conc  = in(:,16);
end


%%% =======================================================================
%%% = END
%%% =======================================================================
