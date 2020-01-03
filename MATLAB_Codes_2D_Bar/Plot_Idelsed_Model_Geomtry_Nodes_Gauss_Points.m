  %%%%%%%%%%%%%%%%%%%%%%%%%%%%  MODULE 02   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%  MODULE 02   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%  MODULE 02   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%  MODULE 02   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
  function Plot_Idelsed_Model_Geomtry_Nodes_Gauss_Points(iGPt,DAT,MDL) 
         Plot_Model_Geometry_wit_Node_Numbers(DAT,MDL)
         Plot_Model_Geometry_with_Cell_Gauss_Points(DAT,MDL) 
         Plot_Active_Influence_Nodes_of_Gauss_Point( iGPt,DAT,MDL)
         Plot_Model_Geometry_with_Boundry_Gauss_Points(DAT,MDL) 
%%------------------------------------------------------- 
%%       Plot the model geomeetry with Node numbers
%%--------------------------------------------------------               
      function Plot_Model_Geometry_wit_Node_Numbers(DAT,MDL) 
  %%  Beam model with its nodes numbered 
           FntS=8;    FntE=8;  FntF=8; FntN=8; NN=25;
           dxN= DAT.Leng/NN;  dyN= DAT.Depth/NN;
           figure(1)
           BeamX=[     0         DAT.Leng     DAT.Leng         0                0      ];
           BeamY=[-DAT.Depth/2  -DAT.Depth/2  DAT.Depth/2  DAT.Depth/2   -DAT.Depth/2  ];
           plot(BeamX, BeamY); hold on; 
         %%  Beam Nodes with its nodes  
           plot( MDL.Node_XCord, MDL.Node_YCord,'o'); hold on;
           for ii=1: MDL.nNodes ;  
               Xp=MDL.Node_XCord(ii)-dxN ; 
               Yp=MDL.Node_YCord(ii)-dyN ;  
               ch1=int2str(ii);  ch=  ['N', ch1];   
               text( Xp,Yp,ch,'color', 'k','FontSize',FntN,'FontWeight','Bold' ); 
           end 
          axis([-.2*DAT.Leng 1.2*DAT.Leng -1.2*DAT.Depth/2 1.2*DAT.Depth/2 ])   
          title(  ' Nodes of Idealised Model')
%%------------------------------------------------------- 
%%       Plot the model geomeetry with Cell and Gauss Points
%%--------------------------------------------------------               
      function Plot_Model_Geometry_with_Cell_Gauss_Points(DAT,MDL) 
           NN=25;  FntN=8;
           figure(2)
           dxG= DAT.Leng/NN;  dyG= DAT.Depth/NN;
           BeamX=[     0         DAT.Leng     DAT.Leng         0                0      ];
           BeamY=[-DAT.Depth/2  -DAT.Depth/2  DAT.Depth/2  DAT.Depth/2   -DAT.Depth/2  ];
           plot(BeamX, BeamY); hold on; 
           for ii=1: MDL.nCells ; 
               xcell(1:4,1)=MDL.Cell_XCord(ii,1:4);
               ycell(1:4,1)=MDL.Cell_YCord(ii,1:4);
               xcell(5,1)=xcell(1,1);ycell(5,1)=ycell(1,1);
               plot(xcell,ycell); hold on;
               xc=sum(xcell(1:4,1))/4; yc=sum(ycell(1:4,1))/4 ;
               ch1=int2str(ii);  ch=  ['C', ch1] ;
               text( xc ,yc ,ch,'color', 'k','FontSize',FntN,'FontWeight','Bold' );  
           end    
           for ii=1:MDL.nGPts
               xp=MDL.XGPt(ii) ;   yp=MDL.YGPt(ii) ; 
               ch1=int2str(ii);  ch=  ['G', ch1]; 
               plot(xp,yp,'X') ;  
               text( xp-dxG ,yp+dyG,ch,'color', 'k','FontSize',FntN,'FontWeight','Bo ld' ); hold on;
           end  
           nn=20; dxN= DAT.Leng/NN;  dyN= DAT.Depth/NN;
           plot( MDL.Node_XCord, MDL.Node_YCord,'o'); hold on;
           for ii=1: MDL.nNodes ;  
               Xp=MDL.Node_XCord(ii)-dxN ; 
               Yp=MDL.Node_YCord(ii)-dyN ;  
               ch1=int2str(ii);  ch=  ['N', ch1];   
               text( Xp,Yp,ch,'color', 'k','FontSize',FntN,'FontWeight','Bold' ); 
           end 
          axis([-.2*DAT.Leng 1.2*DAT.Leng -1.2*DAT.Depth/2 1.2*DAT.Depth/2 ])   
          title(  ' Gauss Points in the cells of Model')    
%%------------------------------------------------------- 
%%       Plot the Active_Inflencing_Nodes_of_Gauss_Point 
%%--------------------------------------------------------               
      function  Plot_Active_Influence_Nodes_of_Gauss_Point( iGPt,DAT,MDL)           
           FntS=8;    FntE=8;  FntF=8; FntN=8; NN=25; NNN=25;
           dxN= DAT.Leng/NN;  dyN= DAT.Depth/NN;
           dxG= DAT.Leng/NNN; dyG= DAT.Depth/NNN;
           figure(3)
           BeamX=[     0         DAT.Leng     DAT.Leng         0                0      ];
           BeamY=[-DAT.Depth/2  -DAT.Depth/2  DAT.Depth/2  DAT.Depth/2   -DAT.Depth/2  ];
           plot(BeamX, BeamY); hold on; 
   %%  Beam Nodes with its nodes  
           plot( MDL.Node_XCord, MDL.Node_YCord,'o'); hold on;
           for ii=1: MDL.nNodes ;  
               Xp=MDL.Node_XCord(ii)-dxN ; 
               Yp=MDL.Node_YCord(ii)-dyN ;  
               ch1=int2str(ii);  ch=  ['N', ch1];   
               text( Xp,Yp,ch,'color', 'k','FontSize',FntN,'FontWeight','Bold' ); 
           end 
         
           del_alp=2*pi/100 ; 
           GXPt=MDL.XGPt(iGPt);GYPt=MDL.YGPt(iGPt);
           plot(GXPt,GYPt,'Xk'); hold on;
           ch1=int2str( iGPt);  ch=  ['G', ch1];
           text( GXPt-dxG,GYPt-dyG,ch,'color', 'k','FontSize',FntN,'FontWeight','Bold' ); 
           for ii=1:100
               alp=ii*del_alp;
               xr(ii) = GXPt + MDL.Serch_Radius*sin(alp);
               yr(ii) = GYPt + MDL.Serch_Radius*cos(alp);
           end
          plot(xr,yr); hold on;
          for ii=1:MDL.GPts_Active_nNodes(iGPt)
              nd=MDL.GPts_Active_Node(iGPt,ii) ; 
              Xp=MDL.Node_XCord(nd)  ; 
              Yp=MDL.Node_YCord(nd)  ;  
              ch1=int2str(nd);  ch=  ['N', ch1]; 
              plot(Xp,Yp,'or') ;  
              text( Xp-dxN,Yp-dyN,ch,'color', 'r','FontSize',FntN,'FontWeight','Bold' ); hold on;
          end 
          axis([-.2*DAT.Leng 1.2*DAT.Leng -1.2*DAT.Depth/2 1.2*DAT.Depth/2 ]) 
          title(  ' Active Nodes of a Gauss Point for Shape Interpolation')
%%------------------------------------------------------- 
%%       Plot the model geomeetry with Boundry Gauss Points
%%--------------------------------------------------------               
      function Plot_Model_Geometry_with_Boundry_Gauss_Points(DAT,MDL) 
           NN=25;  FntN=8;
           figure(4)
           dxG= DAT.Leng/NN;  dyG= DAT.Depth/NN;
           BeamX=[     0         DAT.Leng     DAT.Leng         0                0      ];
           BeamY=[-DAT.Depth/2  -DAT.Depth/2  DAT.Depth/2  DAT.Depth/2   -DAT.Depth/2  ];
           plot(BeamX, BeamY); hold on; 
           for ii=1: DAT.BC_nSides ; 
               for jj=1: MDL.BC_Side_nNodes(ii); 
                   nn=MDL.BC_Side_Node(ii,jj) ;
                   xp=MDL.Node_XCord(nn); yp=MDL.Node_YCord(nn); 
                   ch1=int2str(nn);  ch=  ['N', ch1]; 
                   plot(xp,yp,'o') ;  hold on;
                   text( xp-dxG,yp-dyG,ch,'color', 'k','FontSize',FntN,...
                              'FontWeight','Bold' ); hold on;
               end     
               for jj=1: MDL.BC_Side_nGPts(ii); 
                   nn=MDL.BC_Side_GPt(ii,jj) ;
                   xp=MDL.BC_XGPt(nn); yp=MDL.BC_YGPt(nn); 
                   ch1=int2str(nn);  ch=  ['G', ch1]; 
                   plot(xp,yp,'x') ;  hold on;
                   text( xp-dxG,yp-dyG,ch,'color', 'r','FontSize',FntN,...
                              'FontWeight','Bold' ); hold on;
               end          
           end      
          axis([-.2*DAT.Leng 1.2*DAT.Leng -1.2*DAT.Depth/2 1.2*DAT.Depth/2 ])            
          title(  ' Gauss Points & Nodes on Model Boundry ')