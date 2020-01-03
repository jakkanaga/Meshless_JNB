    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%    Module 03    %%%%%%%%%%%%%%%%%%%%%%% 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%    Module 03    %%%%%%%%%%%%%%%%%%%%%%% 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%    Module 03    %%%%%%%%%%%%%%%%%%%%%%%  
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%    Module 03    %%%%%%%%%%%%%%%%%%%%%%% 
     function  [  PhiMtx ]= Compute_Shape_IntrPolatn_Mtx_Gauss_Point ...
                     (ncell, ngpt, DAT, MDL)
           fprintf('\n  ****  [A],[B],[Phi] of Cell%3i & Gauss Point(local Numb)%3i',...
                 ncell, ngpt) 
           [ nActNds  ActNd  nTerms Serch_Radius  nGPt  xyWJ_GuPt  xyNdPt ]=  ...
                    Get_Cell_GusPoint_Active_Nodes_Loctn(ncell,ngpt,DAT,MDL)
           [ PhiMtx]= Compute_Amtx_Bmtx_Phi_of_Gauss_Point ...
                     (DAT, nGPt,nActNds, ActNd,xyWJ_GuPt,xyNdPt,nTerms,Serch_Radius  ); 
                 
                 
%%--------------------------------------------------------  
%%       Compute_Shape_Functn_of_GausPt_Inside_Boudry
%%--------------------------------------------------------              
     function  [ nActNds  ActNd  nTerms Serch_Radius  nGPt  xyWJ_GuPt  xyNdPt ]=  ...
                    Get_Cell_GusPoint_Active_Nodes_Loctn(ncell,ngpt,DAT,MDL)
            nTerms= DAT.ShFn_nTerms ;   
            nGPt=MDL.Cell_GPt(ncell,ngpt)  
            nActNds =MDL.GPts_Active_nNodes(nGPt)  
            ActNd(1:nActNds,1)= MDL.GPts_Active_Node(nGPt,1:nActNds)   
            xyWJ_GuPt(1)= MDL.XGPt(nGPt) 
            xyWJ_GuPt(2)= MDL.YGPt(nGPt)  
            xyWJ_GuPt(3)= MDL.WGPt(nGPt) 
            xyWJ_GuPt(4)= MDL.JGPt(nGPt)   
            for ii=1: nActNds ;   node=ActNd(ii) ; 
                xyNdPt(ii,1)=MDL.Node_XCord(node);  
                xyNdPt(ii,2)=MDL.Node_YCord(node); 
            end
            Serch_Radius =MDL.Serch_Radius ;
            fprintf('\n  **** Cell_GusPoint_Active_Nodes_Loctn')
            fprintf('\n        Cell: %3i; Gpt %3i; nGPt %3i; xGpt %8.3f; yGpt%8.3f; WGpt %8.3f ; JGpt%8.3f;',...
                 ncell,nGPt, xyWJ_GuPt(1:4))
            fprintf('\n       xNdPt: %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f',...
                  xyNdPt(1:nActNds,1))
            fprintf('\n       yNdPt: %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f',...
                  xyNdPt(1: nActNds,2))  
 
 %%------------------------------------------------------- 
 %%       Compute_AmtxNd_BmtxNd_of_Gauss_Point
 %%------------------------------------------------------- 
     function  [ PhiMtx]= Compute_Amtx_Bmtx_Phi_of_Gauss_Point ...
                     (DAT,nGPt,nActNds, ActNd,xyWJ_GuPt,xyNdPt,nTerms,Serch_Radius  ); 
        fprintf('\n  **** Gauss POint Mtrices Amtx, Bmtx Phi ****')       
        Amtx(1:nTerms,1:nTerms)=0;  dR_dx=1;
        for ii=1:nActNds ;   
            xdif=xyWJ_GuPt(1)-xyNdPt(ii,1) ;
            ydif=xyWJ_GuPt(2)-xyNdPt(ii,2) ;
            Rad_Parm=sqrt(xdif^2+ydif^2)/Serch_Radius ; 
            [ WghtNd(ii) derWght] = Compute_Node_Shape_Weight_Parm(Rad_Parm, dR_dx); 
            PvecNd=[ 1 ; xyNdPt(ii,1) ; xyNdPt(ii,2) ] ; 
            AmtxNd=WghtNd(ii)*PvecNd*PvecNd'  
            Amtx=Amtx+AmtxNd ; 
            BmtxNd(1:nTerms,1)=xyWJ_GuPt(3)*PvecNd(1:nTerms,1) ;
            Bmtx(1:nTerms,ii)=BmtxNd(1:nTerms,1) ; 
        end 
        fprintf('\n Shape Function Weights attached with nodes')
            fprintf('\n   SWght:%8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f',...
                         WghtNd(1:nActNds))
        for ii=1:nTerms
            fprintf('\n    Amtx:%8.3f %8.3f %8.3f %8.3f %8.3f ',Amtx(ii,1:nTerms))
        end
        for ii=1:nTerms
            fprintf('\n    Bmtx:%8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f ',...
                       Bmtx (ii,1:nActNds)) 
        end           
        AInvB_mtx= inv(Amtx)*Bmtx;
        PvecGPt = [ 1 ; xyWJ_GuPt(1) ; xyWJ_GuPt(2)]; 
        Phi=PvecGPt'*AInvB_mtx ; 
        mtx_size=size(Phi); nulmtx=zeros(mtx_size);
        ncol=nActNds*DAT.Each_Node_nUdofs;
        PhiMtx=[ Phi  nulmtx; nulmtx Phi]
        fprintf('\n ***** Phi matrix for Displ at GausPt: %4i ', nGPt) 
        for ii=1:DAT.Each_Node_nUdofs;
            fprintf('\n  Phimtx:%8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f ',...
                       PhiMtx(ii,1:ncol)) 
        end       
%%------------------------------------------------------- 
%%       Compute_Shape_Functn_of_GausPt_Inside_Boudry
%%--------------------------------------------------------             
        function [ Wght derWght] = Compute_Node_Shape_Weight_Parm(Rad_Parm,dR_dx)
            if(Rad_Parm <= 0.5)
                Wght= (2/3)-4.0*Rad_Parm^2+4.0*Rad_Parm^3;
                derWght= (-8.0*Rad_Parm+12.0*Rad_Parm^2)*dR_dx;
            elseif(Rad_Parm <= 1.0)
                Wght= (4/3)-4*Rad_Parm+4*Rad_Parm^2-(4/3)*Rad_Parm^3;
                derWght= (-4+8*Rad_Parm-4*Rad_Parm^2);
            elseif(Rad_Parm > 1.0)
                 Wght=0;
                 derWght=0;
            end
            fprintf('\n  Node_Shape_Weight: %8.3f  %8.3f', Wght,derWght')
%%------------------------------------------------------- 
%%       Compute_Shape_Functn_of_GausPt_Inside_Boudry
%%--------------------------------------------------------                   
            
            
            