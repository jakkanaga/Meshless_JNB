
%%%%%%%%%%%%%%%%%%%%%%%%%%%%  MODULE 00   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%  MODULE 00   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%  MODULE 00   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%  MODULE 00   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%----------------------------------------
%   Main Program of Meshless Beam Model
%---------------------------------------
   function Meshless_Beam_Analysis  
        clear all; clc
        DAT.Title='Meshless_Beam_Analysis  ';
        [DAT  MDL ] = Create_Model_Data_Meshless_Beam_Model( DAT );
        Plot_Idelsed_Model_Geomtry_Nodes_Gauss_Points(12,DAT,MDL) 
          DAT,  MDL    
        ncell=3;   ngpt=4;
      [  PhiMtx ]= Compute_Shape_IntrPolatn_Mtx_Gauss_Point ...
                     (ncell, ngpt, DAT, MDL)
       