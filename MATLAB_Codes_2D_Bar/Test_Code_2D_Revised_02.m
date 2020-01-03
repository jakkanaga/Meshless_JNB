%%%%      2D MATLAB EFG CODE
%%------------------------------------------------------- 
%%        Test_Code_2D_Revised
%%--------------------------------------------------------
     function Test_Code_2D_Revised
          clear all; clc
%%   Data structure :  DAT   Basic Data of the Simulation Model  
%%   Data structure :  MDL   Model Data Created of Simulation Model  
           [ DAT ]   = Basic_Model_Data_Description  
         [DAT  MDL ] = Create_Nodal_Data_Mesh_Grid_Model( DAT );
         [DAT  MDL ] = Select_Gauss_Qudrature_Two_Points ( DAT, MDL ) ;
         [DAT  MDL ] = Create_Cell_Mesh_Grid_Points_Model (DAT, MDL ) ;
         [DAT  MDL ] = Compute_Gauss_Points_with_Parm_Model(DAT, MDL)
         [DAT  MDL ] = Compute_Active_Nodes_for_Gauss_Quad_Points(DAT,MDL); 
         [DAT  MDL ] = Compute_Active_Nodes_for_Gauss_Quad_Points(DAT,MDL); 
         [DAT  MDL ] = Compute_Nodes_on_Displ_Boundary_Model(DAT,MDL) 
 
%%------------------------------------------------------- 
%%           Basic_Model_Data_Description  
%%--------------------------------------------------------
     function [ DAT ]= Basic_Model_Data_Description 
          DAT.Leng = 2;                               %  Beam Length
          DAT.Depth=2;                                %  Beam Depth
          DAT.YMod = 30e6;                            %  Youngs Modulus
          DAT.Nu=0.3;                                 %  Poission Ratio
          DAT.Tip_Load=1000;                          %  Tip Load
          DAT.nCell_Leng =4;                          %  # cells on Length
          DAT.nCell_Depth = 1;                        %  # cells on Depth         
          DAT.Serch_NDP=1.5;                          %  NonDim Search_Parm
          DAT.GSchm_nPts=2                            %  Gauss Scheme-#Pts 
        
           
%%------------------------------------------------------- 
%%       Create_Nodal_Data_Mesh_Grid_Model
%%--------------------------------------------------------                  
    function [DAT  MDL ] = Create_Nodal_Data_Mesh_Grid_Model( DAT )           
%% SET UP NODAL COORDINATES
          MDL.nNodes=(DAT.nCell_Leng+1)*(DAT.nCell_Depth+1);%  # nodes Model
          MDL.Cell_Depth = DAT.Depth /DAT.nCell_Depth;      %  Cell depth  
          MDL.Cell_Leng = DAT.Leng /DAT.nCell_Leng;         %  Cell Length
          for ii = 1:DAT.nCell_Leng+1  
              for jj = 1 : DAT.nCell_Depth+1;
                  nd= (ii-1)*(DAT.nCell_Depth+1)+jj;
     %% Coordinates of the node             
                  MDL.Node_XCord(nd,1) = (ii-1)*MDL.Cell_Leng ;
                  MDL.Node_YCord(nd,1) = DAT.Depth/2 -(jj-1)*MDL.Cell_Depth ;
               end
          end   
          fprintf('\n  Nodal XCord:  %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f', ...
                 MDL.Node_XCord(1:MDL.nNodes,1))  
          fprintf('\n  Nodal YCord:  %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f', ...
                 MDL.Node_YCord(1:MDL.nNodes,1))    
          
%%------------------------------------------------------- 
%%         Create_Cell_Mesh_Grid_Points_Model 
%%--------------------------------------------------------                  
    function [DAT  MDL ] =Create_Cell_Mesh_Grid_Points_Model (DAT, MDL ) 
         
          MDL.nCells=DAT.nCell_Leng*DAT.nCell_Depth;        %  # Cells Model      
          MDL.Cell_nGQPts=DAT.GSchm_nPts^2;                  %  # cell GQadPts
          for ii = 1:DAT.nCell_Leng   
              for jj = 1 : DAT.nCell_Depth 
                  nd(1)= (ii-1)*(DAT.nCell_Depth+1)+jj;
                  nd(2)= nd(1)+(DAT.nCell_Depth+1);
                  nd(3)=nd(2)+1  ;  nd(4)=nd(1)+1;
                  ncell=(ii-1)*DAT.nCell_Depth +jj;
                  MDL.Cell_Nodes(ncell,1:4)=nd(1:4);  % Nodesl of the cells
                  for kk=1:MDL.Cell_nGQPts  
    %%  Cooedinates of the of each cell of the model                    
                      MDL.Cell_XCord(ncell,kk)=MDL.Node_XCord(nd(kk));
                      MDL.Cell_YCord(ncell,kk)=MDL.Node_YCord(nd(kk));
                  end                       
               end
          end 
          fprintf('\n  Create_Cell_Mesh_Grid_Points_Model')  
          for ii=1: MDL.nCells 
              fprintf('\n  %3i- th Cell  Nodes:%8i %8i %8i %8i', ...
                        ii, MDL.Cell_Nodes(ii,1:MDL.Cell_nGQPts))
              fprintf('\n  %3i- th Cell XCord:  %8.3f %8.3f %8.3f %8.3f', ...
                        ii, MDL.Cell_XCord(ii,1:MDL.Cell_nGQPts))
              fprintf('\n  %3i -th Cell YCord:  %8.3f %8.3f %8.3f %8.3f', ...
                        ii, MDL.Cell_YCord(ii,1:MDL.Cell_nGQPts))  
          end
          
%%------------------------------------------------------- 
%%         Select_Gauss_Qudrature_Four_Points
%%--------------------------------------------------------                  
        function  [DAT  MDL ] = Select_Gauss_Qudrature_Four_Points ( DAT, MDL )      
%%  Gauss quadrature: 4 points and its weights  for limits of(-1 to1)          
           DAT.GCPt(1,1) = -0.861136311594052575224;
           DAT.GCPt(2,1) = -0.339981043584856264803;
           DAT.GCPt(3,1) = -DAT.GCPt(2,1);
           DAT.GCPt(4,1) = -DAT.GCPt(1,1);
           DAT.GCWt(1,1) = 0.347854845137453857373;
           DAT.GCWt(2,1) = 0.652145154862546142627;
           DAT.GCWt(3,1) = DAT.GCWt(1,1);
           DAT.GCWt(4,1) = DAT.GCWt(2,1);
           fprintf('\n  GPt:%8.3f %8.3f %8.3f %8.3f',DAT.GCPt(1:4,1))             
           fprintf('\n  Gwt:%8.3f %8.3f %8.3f %8.3f',DAT.GCWt(1:4,1))  
%%------------------------------------------------------- 
%%          Select_Gauss_Qudrature_Two_Points
%%--------------------------------------------------------                  
        function  [DAT  MDL ] = Select_Gauss_Qudrature_Two_Points ( DAT, MDL )      
%%  Gauss quadrature: 4 points and its weights  for limits of(-1 to1)          
           DAT.GCPt(1,1) =  -0.577350269
           DAT.GCPt(2,1) =   0.577350269 
           DAT.GCWt(1,1) =   1.0
           DAT.GCWt(2,1) =   1.0
           fprintf('\n  GPt:%8.3f %8.3f %8.3f %8.3f',DAT.GCPt(1:2,1))             
           fprintf('\n  Gwt:%8.3f %8.3f %8.3f %8.3f',DAT.GCWt(1:2,1))       
%%------------------------------------------------------- 
%%      Compute_Gauss_Points_with_Parm_Model
%%--------------------------------------------------------                 
     function [ DAT, MDL ] = Compute_Gauss_Points_with_Parm_Model(DAT, MDL)
 %%    number of gauss points in the cell 
        nGPts=0  
        for ic=1:MDL.nCells
 %%    coordinates of the cell corners           
           Cell_Xcor(1:4,1) = MDL.Cell_XCord(ic,1:MDL.Cell_nGQPts);   % Xcor of cell vertices
           Cell_Ycor(1:4,1) = MDL.Cell_YCord(ic,1:MDL.Cell_nGQPts);  % Ycor of cell vertices
  %%   Loops for  loop for Gpt in X-dir & Y-dir       
           for jc=1:DAT.GSchm_nPts ;              
               for kc=1:DAT.GSchm_nPts 
                 nGPts = nGPts+1;
  %%   psi and eta isoparametric coordinates and weights                
                     p = DAT.GCPt(jc,1) ;  
                     e = DAT.GCPt(kc,1) ;  % QPt Gauss pt                
                    wp = DAT.GCWt(jc,1) ; 
                    we = DAT.GCWt(kc,1) ;   % WPt Gauss Wt
  %%   isoparametric coordinates to quadrilateral coordinates transfmtn
                    TMtx = 0.25*[(1-p)*(1-e) (1+p)*(1-e) (1+p)*(1+e) (1-p)*(1+e)] ; 
  %%   Gauss points in the quadrilateral and its weight  
                    MDL.XCPt(nGPts) =TMtx*Cell_Xcor ; 
                    MDL.YCPt(nGPts) =TMtx*Cell_Ycor ; 
                    MDL.WCPt(nGPts) = wp*we  ;
  %%   derevetives  of coordinate transformation  and Jacobian                
                    dTbydp =0.25*[ -(1-e) +(1-e) +(1+e) -(1+e) ] ;
                    dTbyde =0.25*[ -(1-p) -(1+p) +(1+p) +(1-p) ] ; 
                    dxbydp = dTbydp*Cell_Xcor  ;
                    dxbyde = dTbyde*Cell_Xcor  ;
                    dybydp = dTbydp*Cell_Ycor  ;
                    dybyde = dTbyde*Cell_Ycor  ;
                    MDL.CJacbn(nGPts) = dxbydp*dybyde -dxbyde*dybyde ;
               end
           end 
        end    
           fprintf('\n   Parameters of Guassian Quadrature Points') 
           fprintf('\n     #       Xcor     Ycor     Weight     Jacobian')
           for ii=1:nGPts
               fprintf('\n  %4i   %8.3f  %8.3f  %8.3f  %8.3f',ii, MDL.XCPt(ii),...
                         MDL.YCPt(ii),MDL.WCPt(ii),MDL.CJacbn(ii))
           end 
           MDL.nTGPts=nGPts ;  %   Total Gauss points of model
%%------------------------------------------------------- 
%%       Compute Nodes Inflencing Gauss Quad Points
%%--------------------------------------------------------              
     function [ DAT MDL ]= Compute_Active_Nodes_for_Gauss_Quad_Points(DAT,MDL) 
        MDL.Serch_Radius=MDL.Cell_Leng*DAT.Serch_NDP; %  Search Radius
        if(MDL.Cell_Leng > MDL.Cell_Depth)  
            MDL.Serch_Radius= MDL.Cell_Depth*DAT.Serch_NDP;
        end    
        for ii=1: MDL.nTGPts
            XGPt= MDL.XCPt(ii); YGPt= MDL.YCPt(ii); 
            XDiff = abs(MDL.Node_XCord - XGPt*ones(MDL.nNodes,1)) ;
            YDiff = abs(MDL.Node_YCord - YGPt*ones(MDL.nNodes,1)) ;
            Radius= sqrt(XDiff.*XDiff+YDiff.*YDiff);
            Radius_Diff = abs(MDL.Serch_Radius- Radius) ;
            Active_Node= find(Radius_Diff < MDL.Serch_Radius); 
            nn=size(Active_Node,1);
            MDL.GPts_Active_nNodes(ii)=nn;
            MDL.GPts_Active_Node(ii,1:nn)=Active_Node(1:nn,1);
        end
        fprintf('\n   Gauss Quad Points and its Active Nodes')
        fprintf('\n   Gpt  #AcNds     Active Nodes')
        for ii=1: MDL.nTGPts
            nn=MDL.GPts_Active_nNodes(ii); 
            Active_Node(1:nn,1)=MDL.GPts_Active_Node(ii,1:nn); 
            fprintf('\n  %4i %4i    %4i %4i %4i %4i %4i %4i %4i %4i %4i %4i %4i %4i', ...
                       ii,nn, Active_Node(1:nn,1) )
       end  
           
%%------------------------------------------------------- 
%%     Compute_Nodes_on_Displ_Boundary_Model
%%--------------------------------------------------------              
     function [ DAT MDL ]= Compute_Nodes_on_Displ_Boundary_Model(DAT,MDL) 
  %%  Nodes on the top of the beam surface
        nTop_Nodes=0 ; nBot_Nodes=0 ; 
        for ii=1:MDL.nNodes;
            if (MDL. Node_YCord(ii)==DAT.Depth/2);
                nTop_Nodes=nTop_Nodes+1;
                Top_Node(nTop_Nodes)=ii ;   
            end  
            if (MDL. Node_YCord(ii)==-DAT.Depth/2);
                nBot_Nodes=nBot_Nodes+1;
                Bot_Node(nBot_Nodes)=ii ;   
            end
        end
        MDL.BC_nTop_Nodes = nTop_Nodes;
        MDL.BC.Top_Node=Top_Node(1: nTop_Nodes);
        MDL.BC_nBot_Nodes = nBot_Nodes;
        MDL.BC.Bot_Node=Bot_Node(1: nBot_Nodes);
        fprintf('\n   Displ-Boudry Nodes lying on Beams Surface ')
        fprintf('\n   Top Beam Surface: # Number of Nodes: %4i and their Details',...
                        MDL.BC_nTop_Nodes)
        fprintf('\n  %4i %4i %4i %4i %4i %4i %4i %4i %4i %4i %4i %4i %4i %4i', ...
                       MDL.BC.Top_Node(1: MDL.BC_nTop_Nodes))            
        fprintf('\n   Bot Beam Surface: # Number of Nodes: %4i and their Details',...
                        MDL.BC_nBot_Nodes)
        fprintf('\n  %4i %4i %4i %4i %4i %4i %4i %4i %4i %4i %4i %4i %4i %4i', ...
                       MDL.BC.Bot_Node(1: MDL.BC_nBot_Nodes))    
           
%%------------------------------------------------------- 
%%       Compute_Gaus_Points_on_Tractn_Boudry_Model
%%--------------------------------------------------------              
     function [ DAT MDL ]= Compute_Gaus_Points_on_Tractn_Boudry_Model(DAT,MDL) 
  %%  Nodes on the top of the beam surface
        nTop_Nodes=0 ; nBot_Nodes=0 ; nRSide_Nodes=0
        for ii=1:MDL.nNodes;
            ii, Node_XCord(ii)
            if (MDL. Node_YCord(ii)==DAT.Depth/2);
                nTop_Nodes=nTop_Nodes+1;
                Top_Node(nTop_Nodes)=ii ;   
            end  
            if (MDL. Node_YCord(ii)==-DAT.Depth/2);
                nBot_Nodes=nBot_Nodes+1 
                Bot_Node(nBot_Nodes)=ii ;   
            end
%             MDL.Node_XCord
%             if (MDL.Node_XCord(ii)==DAT.Leng);
%                 nRSide_Nodes=nRSide_Nodes+1 
%                 nRSide_Node(nRSide_Nodes)=ii    
%             end
%             if (MDL.Node_XCord(ii)==0.);
%                 nLSide_Nodes=LSide_Nodes+1 
%                 nLSide_Node(LSide_Nodes)=ii    
%             end
        end
        MDL.BC_nTop_Nodes = nTop_Nodes;
        MDL.BC_Top_Node=Top_Node(1: nTop_Nodes);
        MDL.BC_nBot_Nodes = nBot_Nodes;
        MDL.BC_Bot_Node=Bot_Node(1: nBot_Nodes);
        MDL.BC_RSide_nNodes = nRSide_Nodes 
        MDL.BC_nRSide_Node=nRSide_Node(1: nLSide_Nodes) 
        MDL.BC_nLSide_Nodes = nLSide_Nodes 
        MDL.BC_nLSide_Node=nLSide_Node(1: nLSide_Nodes) 
        fprintf('\n   Boundary Nodes lying on Beams Surfaces ')
        fprintf('\n   Top Beam Surface: # Number of Nodes: %4i and their Details',...
                        MDL.BC_nTop_Nodes)
        fprintf('\n  %4i %4i %4i %4i %4i %4i %4i %4i %4i %4i %4i %4i %4i %4i', ...
                       MDL.BC_Top_Node(1: MDL.BC_nTop_Nodes))            
        fprintf('\n   Bot Beam Surface: # Number of Nodes: %4i and their Details',...
                        MDL.BC_nBot_Nodes)
        fprintf('\n  %4i %4i %4i %4i %4i %4i %4i %4i %4i %4i %4i %4i %4i %4i', ...
                       MDL.BC_Bot_Node(1: MDL.BC_nBot_Nodes))             
        fprintf('\n   RSide Beam Surface: # Number of Nodes: %4i and their Details',...
                        MDL.BC_RSide_nNodes)
        fprintf('\n  %4i %4i %4i %4i %4i %4i %4i %4i %4i %4i %4i %4i %4i %4i', ...
                       MDL.BC_RSide_Node(1: MDL_BC_RSide_Nodes))      
        fprintf('\n   LSide Beam Surface: # Number of Nodes: %4i and their Details',...
                        MDL.BC_LSide_Node)
        fprintf('\n  %4i %4i %4i %4i %4i %4i %4i %4i %4i %4i %4i %4i %4i %4i', ...
                       MDL.BC_LSide_Node(1: MDL_BC_LSide_Nodes))       
           
           
           
           
           
                       
% % SET UP CONNECTIVITY ARRAY
%           for j =1 :ndivl
%               for i = 1:ndivw
%                   elemn = (j-1)*ndivw + i;
%                   nodet(elemn,1) = elemn + (j-1) ;
%                   nodet(elemn,2) = nodet (elemn, 1) + 1;
%                   nodet(elemn,3) = nodet (elemn, 2) +ndivw+1 ;
%                   nodet(elemn,4) = nodet (elemn,3) -1 ;
%               end
%              
%              
%              
%              
%              
%              
%              
%              
%              
%              
%              
%              
%           k = sparse (numnod*2 ,numnod*2) ;
%           for gg=gs
%               gpos = gg(1:2);
%               weight = gg(3);
%               jac = gg(4);        
% % DETERMINE NODES IN NEIGHBORHOOD OF GAUSS POINT
% 
%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  PAGE 235  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%               v = domain(gpos ,x ,dm,numnod) ;
%               L = length(v);
%               en = zeros (1, 2*L) ;
%               [phi, dphix, dphiy] = shape (gpos, dmax, x, v , dm) ;
%               Bmat=zeros(3 ,2*L);
%               for j=1:L
%                   Bmat(1:3,(2*j-1):2*j) = ...
%                       [dphix(j) 0;0 dphiy(j) ;dphiy(j) dphix(j)];
%               end
%               for i=1:L
%                   en(2+i-1) = 2*v(i)-1;
%                   en(2*i) = 2*v(i) ;
%               end
%          
%               k(en,en) = k(en,en)+sparse((weight*jac)*Bmat'*Dmat*Bmat) ;
%            end
% % DETERMINE NODES ON BOUNDARY, SET UP BC'S
%            ind1 = 0;ind2 = 0;
%            for j=l:numnod
%                if (x(1, j)==0.0)
%                    indl=indl+l;
%                    nnu(1,ind1) = x(1,j);
%                    nnu(2,ind1) = x(2,j);
%                end
%                if (x(1, j)==Lb)
%                    ind2=ind2+1;
%                    nt(1,ind2) = x(1,j);
%                    nt(2,ind2) = x(2,j);
%                end
%            end
%            lthu = length(nnu) ;
%            ltht = length(nt);
%            ubar = zeros(lthu*2,1);
%            f = zeros (nwrniod*2,1);
% %SET UP GAUSS POINTS ALONG TRACTION BOUNDARY
%            i1id=0 ;
%            gauss=pgauss(quado);
%            for i=1: (ltht-1)
%                ycen= (nt (2, i+1)+nt(2,i))/2 ;
%                jcob =abs((nt(2,i+l)-nt(2,i)) /2) ;
%                for j=l:quado
%                    mark(j) = ycen-gauss(1,j)*jcob;
%                    ind = ind+1;
%                    gst(1,ind) = nt(1,i);
%                    gst(2,ind) = mark(j) ;
%                    gst(3,ind) = gauss(2,j) ;
%                    gst(4,ind) = jcob;
%                end
%            end
%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  PAGE 236  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%            
% %SET UP GAUSS POINTS ALONG DISPLACEMENT BOUNDARY
%             gsu = gst
%             gsu(1,1:i) = zeros(1,ind);
%             qk = zeros(1:2*lthu);
% %INTEGRATE FORCES ALONG BOUNDARY
%             Imo = (1/12)*Dn3;
%             for gt=gst
%                 gpos = gt(1:2) ;
%                 weight= gt(3) ;
%                 jac = gt(4);
%                 v = domain(gpos,x,dm,numnod) ;
%                 L = length(v) ;
%                 en = zeros(1,2*L);
%                 force=zeros (1,2*L) ;
%                 [phi,dphix,dphiy] = shape(gpos,dmax,x,v,dm);
%                 tx = 0 ;
%                 ty = -(P/(2*Imo))*((D^2)/4-gpos(2,1)^2) ;
%                 for i=l:L
%                      en(2*i-1) = 2*v(i)-1;
%                      en(2*i)=2*v(i) ;
%                      force(2*i-l)=tx*phi(i);
%                      force(2*i) = ty*phi (i) ;
%                 end
%                 f(en) = f(en)+jac*weight*force' ;
%            end
% % INTEGRATE G MATRIX AND Q VECTOR ALONG DISPLACEMENT BOUNDARY
%            GG = zeros(numnod*2,lthu*2);
%            indl = 0;   ind2=0;
%            for i=l:(lthu-l)
%                indl=indl+l;
%                ml = indl; m2 = ml+l;
%                yl = nnu(2,mi) ; y2 = nnu(2,m2) ;
%                len = yl-y2;
%                for j=l:quado
%                    ind2=ind2+1;
%                    gpos = gsu(l:2,ind2) ;
%                    weight = gsu(3,j);
%                    jac = gsu(4,j);
%                    facll = (-P*nnu(2,mi)) / (6*young*Imo) ;
%                    fac2 = P/ (6*young*Imo) ;
%                    xpl = gpos(1,1);
%                    ypl = gpos(2,l) ;
%                    uxexl = (6*Lb-3*xpl)*xp1+(2+nu)*(yp1^2-(D/2)^2) ;
%                    uxexl = uxexl*facll;
%                    uyexl = 3*nu*ypln2*(~b-xpl)+0.25*(4+5*nu)*xpl*~n2+(3*Lb-xpl)*xp1n2;
%                    uyexl = uyexl*fac2;
%  
%   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  PAGE 237  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%                   
%                    
%                    Nl = (gpos(2,I)-y2)/len; N2 = 1-N1;
%                    qk(2*ml-I) = qk(2*ml-I) -weight* jac*Nl*uxexl ;
%                    qk(2*ml) = qk(2*ml) - weight* jac*NI*uyexl ;
%                    qk(2*m2-1) = qk(2*m2-1) - weight*jac*N2+uxexl;
%                    qk(2*m2) = qk(2*m2) - weight*jac*N2*uyexl;
%                    v = domain(gpos ,x ,dm,numnod) ;
%                    [phi,dphix,dphiy] = shape(gpos,dmax,x,v,dm);
%                    L = length(v);
%                    for n=l:L
%                        GI = -weight*jac*phi(n)*[Nl 0;0 N1] ;
%                        G2 = -weight*jac*phi(n)*[N2 0;0 N2];
%                        cl=2*v(n)-l;   c2=2*v(n);   c3=2*m1-l;   c4=2*m1;
%                        c5=2*m2-l;c6=2*m2;
%                        GG(c1:c2,c3:c4) = GG(c1:c2,c3:c4)+G1;
%                        GG(c1:c2,c5:c6) = GG(c1:c2,c5:c6)+ G2 ;
%                    end
%                 end
%            end
%                
% % ENFORCE BC'S USING LAGRANGE MULTIPLIERS
%             f = [f;zeros(lthu*2,1)] ;
%             f(numnod*2+1:numnod*2+2*lthu,l) = -qk';
%             m = sparse([k GG;GG' zeros(lthu*2)]);
%             d = m\f ;
%             u = d(l:2*numnod) ;
%             for i=l:numnod
%                 u2(1,1) = u(2*i-I) ;
%                 u2(2,i) = u(2*i);
%             end
% % SOLVE FOR OUTPUT VARIABLES - DISPLACEMENTS
%             for gg=x
%                 ind = ind+l;
%                 gpos = gg(1:2);
%                 v = domain(gpos,x,dm,numod) ;
%                 [phi,dphix,dphiy] = shape (gpos,dmax,x,v,dm) ;
%                 displ(2*ind-1) = phi*u2(l,v)' ;
%                 displ(2*ind) = phi*u2(2,v)' ;
%             end
% % SOLVE FOR STRESSES AT GAUSS POINTS
%             ind = 0;
%             for gg=gs
%                 ind = ind+l;
%                 gpos = gg(l:2) ;                   
%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  PAGE 238  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%  
%                 weight = gg(3);
%                 jac = gg(4);
%                 v = domain(gpos,x,dm,numnod) ;
%                 L = length(v) ;
%                 en = zeros (I, 2*L) ;
%                 [phi, dphix , dphiy] = shape (gpos , dmax , x , v , dm) ;
%                 Bmat=zeros (3,2*L) ;
%                 for j=l:L
%                     Bmat(l:3, (2*j-I) :2*j) = [dphix(j) 0;O dphiy(j) ;dphiy(j) dphix(j)] ;
%                 end
%                 for i=l:L
%                     en(2*i-I) = 2*v(i)-I;
%                     en(2*i) = 2*v(i);
%                 end
%                 stress(l:3,ind) = Dmat*Bmat*u(en);
%                 stressex(1 ,ind) = (l/Imo)*P*(Lb-gpos(1,I) )*gpos(2,1) ;
%                 stressex (2, ind) = 0 ;
%                 stressex(3, ind) = -0.5* (P/Imo) * (0.25*DA2 - gpos (2, I) -2) ;
%                 err = stress(l:3,ind)-stressex(l:3,ind);
%                 err2 = weight*jac(0.5*(inv(Dmat)*err)'*(err)) ;
%                 enorm = enorm + err2;
%             end
%             enorm = sqrt (enorm)
%             
%             
%  
% % MESH GENERATION PROGRAM FOR 2D BEAM IN BENDING
%        function[x conn numcell numq] = mesh2(length,height,ndivl,ndivw)
% % INPUT DATA
%             numcell= ndivw*ndivl;
%             numq = (ndivl+1) * (ndivw+1) ;
% % SET UP NODAL COORDINATES
%            for i = 1:(ndivl+1)
%                for j = 1 : (ndivw+1)
%                    x(1,((ndivw+1)*(i-1) +j))= (length/ndivl)*(i-1);
%                    x(2,((ndivw+1)*(i-1) +j))= -(height/ndivw)*(j-1)+height/2;
%                end
%            end
% % SET UP CONNECTIVITY ARRAY
%           for j =1 :ndivl
%               for i = 1:ndivw
%                   elemn = (j-1)*ndivw + i;
%                   nodet(elemn,1) = elemn + (j-1) ;
%                   nodet(elemn,2) = nodet (elemn, 1) + 1;
%                   nodet(elemn,3) = nodet (elemn, 2) +ndivw+1 ;
%                   nodet(elemn,4) = nodet (elemn,3) -1 ;
%               end
% %  
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  PAGE 238  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%  
%           end
%             conn = nodet';
%          
% % routine to set up gauss points, jacobian, and weights              
%      function [gs] = egauss(xc,conn,gauss,numcell)
%            index=0 ;
%            one = ones(1,4);
%            psiJ = [-1,+1,+1,-1]; etaJ = [-1,-1,+1,+1];
%            l = size(gauss);
%            l = l(2);
%            for e=l:numcell
% % DETERMINE NODES IN EACH CELL
%               for j = 1:4
%                   je=conn(j,e) ;  xe(j)=xc(1,je) ;  ye(j)=xc(2, je) ;
%               end
%               for i=l:l
%                   for j=i:l
%                       index = index+1;
%                       eta=gauss(1,i);      psi=gauss(1, j) ;
%                       N = .25*(one+psi*psiJ).*(one+eta*etaJ) ;
%                       NJpsi = .25*psiJ.*(one+eta*etaJ) ;
%                       NJeta=.25*etaJ.*(one+psi*psiJ) ;
%                       xpsi=NJpsi*xe';   ypsi=NJpsi*ye';
%                       xeta=NJeta*xe';   yeta=NJeta*ye';
%                       jcob=xpsi*yeta-xeta*ypsi;
%                       xq = N*xe';    yq = N*ye';
%                       gs(1,index) = xq;
%                       gs(2,index) = yq;
%                       gs(3,index) = gauss(2,i)*gauss(2,j) ;
%                       gs(4,index) = jcob;
%                   end
%               end
%            end
%        
%  %%%  -----------------------------------------------------------------
%  % This function returns a matrix with 4 gauss points and their weights
%  %%% ------------------------------------------------------------------
%       function v = pgauss(k)
%            v(1,1) =-.861136311594052575224;
%            v(1,2) =-.339981043584856264803;
%            v(1,3) = -v(1,2);
%            v(1,4) = -v(1, 1) ;
%            v(2,1) =.347854845137453857373;
%            v(2,2) =.652145154862546142627;
%            v(2,3) = v(2,2);
%            v(2,4) = v(2,1);
% %%%%% ---------------------------------------------------
% %  DETERMINES NODES IN DOMAIN OF INFLUENCE OF POINT GPOS
% %%%%% ---------------------------------------------------
%      function v = domain(gpos,x,dm,nnodes)
%           dif = gpos*ones(1,nnodes)-x;
%           a = dm-abs (dif) ;
%           b = [a;zeros(1 ,nnodes)] ;
%           c = all (b>=-100*eps) ;
%           v = find(c) ;
% %%%%%  page 240  %%%%%%
% %%%%% ---------------------------------------------------
% % EFG shape function and it's derivative with linear base
% %%%%% ---------------------------------------------------
%      function [ phi dphix dphiy] = shape(gpos,dmax,x,v,dm)
% % EFG shape function and it's derivative with linear base
%          L = length(v) ;
%          won = ones(1,L);
%          nv = x(1:2,v);
%          p = [won;nv] ;
%          dif = gpos*won-nv;
%          t = dm(1:2,v)/dmax;
% % WEIGHTS--W and dw are vectors
%          [w , dwdx , dwdy] = cubwgt (dif , t , v , dmax ,dm) ;
%          B = p.*[w;w;w];
%          pp = zeros(3) ;
%          aa = zeros (3) ;
%          daax = zeros(3);
%          daay = zeros (3) ;
%          for i=1:L
%              pp = p(1:3,i)*p(1:3,i)' ;
%              aa = aa+w(1,i)*pp;
%              daax = daax+dwdx(1,i)*pp;
%              daay = daay+dwdy(1,i)*pp;
%          end
%         pg =[ 1  gpos'] ;
%         [L,U,PERM] = lu(aa);
%         for i=1:3
%             if i==1
%                C = PERM*pg';
%             elseif i==2
%                C = PERM*([0 1 0]' - daax*gam(1:3,1));
%             elseif i==3
%                C = PERM*([0 0 1]' - daay*gam(1:3,1));
%             end
%             D1 = C(1);
%             D2 = (C(2) - L(2,1)*D1) ;
%             D3 = (C(3) - L(3,1)*D1-L(3,2)*D2) ;
%             gam(3,i) = D3/U(3,3) ;
%             ga1(2,i) = (D2-U(2,3)*gam(3,i))/(U(2,2));
%             gam(1 ,i) =(D1-U(1,2)*gam(2,i)-U(1,3)*gam(3,i))/U(1,1) ;
%         end
%         phi = gam(1:3,1)'*B;
%         dbx = p.* [ dwdx; dwdx; dwdx] ;
%         dby = p.* [ dwdy; dwdy; dwdy] ;
%         dphix = gam(1:3,2)'*B + gam(1:3,1)'*dbx;
%         dphiy = gam(1:3,3)'*B + gam(1:3,1)'*dby;
% %%%%% ---------------------------------------------------       
%            % CUBIC SPLINE WEIGHT FUNCTION  
% %%%%% --------------------------------------------------- 
%     function [w , dwdx , dwdy] = cubwgt (dif , t , v , dmax , dm)
%           l = length(v) ;
%           for i=1:l
% %%%%    page 241   %%%%%%%%%  
%               drdx=sign(dif(1,i))/dm(1,v(i));
%               drdy=sign(dif(2,i))/dm(2,v(i));
%               rx=abs(dif(1,i))/dm(1,v(i));
%               ry=abs(dif(2,i))/dm(2,v(i));
%               if (rx > 0.5)
%                  wx=(4/3)-4*rx+4*rx*rx-(4/3)*rx^3; 
%                  dwx=(-4+8*rx+4*rx^2)*drdx;
%               end
%               if (rx <= 0.5)
%                  wx=(2/3)-4*rx*rx+4*rx^3; 
%                  dwx=(-8*rx+12*rx^2)*drdx;
%               end   
%               if (ry > 0.5)
%                  wy=(4/3)-4*ry+4*ry*ry-(4/3)*ry^3; 
%                  dwy=(-4+8*ry-4*ry^2)*drdy;
%               end
%               if (rx <= 0.5)
%                  wy=(2/3)-4*ry*ry+4*ry^3; 
%                  dwy=(-8*ry+12*ry^2)*drdy;
%               end  
%               w(i)=wx*wy;
%               dwdx(i)=wy*dwx;
%               dwdy(i)=wx*dwy
%           end      
%  