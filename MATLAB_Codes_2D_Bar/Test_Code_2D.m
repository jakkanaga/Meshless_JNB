%%%%      2D MATLAB EFG CODE - SQUARE DOMAIN
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  PAGE 234  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
clear ;
% DEFINE BOUNDARIES/PARAMETERS
          Lb = 48;
          D=12;
          young = 30e6;   nu=0.3;
          P=1000 ;
% PLANE STRESS DMATRIX
          Dmat = (young/(1-nu^2))*[ 1 nu 0; nu 1 0; 0  0 (1-nu)/2] ;
% SET UP NODAL COORDINATES
          ndivl = 10;   ndivw = 4;
         [x,conn1,conn2,numnod] = mesh2 (Lb,D,ndivl,ndivw);
% DETERMINE DOMAINS OF INFLUENCE - UNIFORM NODAL SPACING
          dmax=3.5 ;
          xspac = Lb/ndiv1;
          yspac = D/ndivw;
          dm(1,1:numnod)=dmax*xspac*ones(1,numnod);
          dm(2,1:numnod)=dmax*yspac*ones(1,numnod) ;
% SET UP QUADRATURE CELLS
          ndivlq = 10;ndivwq = 4;
          [xc,conn,numcell,numq] = mesh2(Lb,D,ndivlq,ndivwq);
% SET UP GAUSS POINTS, WEIGHTS, AND JACOBIAN FOR EACH CELL
          quado = 4;
          [gauss] = pgauss (quado) ;
          numq2 = numcell*quad^2;
          gs = zeros (4,numq2) ;
          [gs] = egauss(xc,conn,gauss,numcell);
% LOOP OVER GAUSS POINTS TO ASSEMBLE DISCRETE EQUATIONS
          k = sparse (numnod*2 ,numnod*2) ;
          for gg=gs
              gpos = gg(1:2);
              weight = gg(3);
              jac = gg(4);        
% DETERMINE NODES IN NEIGHBORHOOD OF GAUSS POINT

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  PAGE 235  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

              v = domain(gpos ,x ,dm,numnod) ;
              L = length(v);
              en = zeros (1, 2*L) ;
              [phi, dphix, dphiy] = shape (gpos, dmax, x, v , dm) ;
              Bmat=zeros(3 ,2*L);
              for j=1:L
                  Bmat(1:3,(2+j-1):2*j) = [dphix(j) 0;0 dphiy(j) ;dphiy(j) dphix(j)];
              end
              for i=l:L
                  en(2+i-1) = 2*v(i)-1;
                  en(2*i) = 2*v(i) ;
              end
              k(en,en) = k(en,en)+sparse((weight*jac)*Bmat'*Dmat*Bmat) ;
           end
% DETERMINE NODES ON BOUNDARY, SET UP BC'S
           ind1 = 0;ind2 = 0;
           for j=l:numnod
               if (x(1, j)==0.0)
                   indl=indl+l;
                   nnu(1,ind1) = x(1,j);
                   nnu(2,ind1) = x(2,j);
               end
               if (x(1, j)==Lb)
                   ind2=ind2+1;
                   nt(1,ind2) = x(1,j);
                   nt(2,ind2) = x(2,j);
               end
           end
           lthu = length(nnu) ;
           ltht = length(nt);
           ubar = zeros(lthu*2,1);
           f = zeros (nwrniod*2,1);
%SET UP GAUSS POINTS ALONG TRACTION BOUNDARY
           i1id=0 ;
           gauss=pgauss(quado);
           for i=1: (ltht-1)
               ycen= (nt (2, i+1)+nt(2,i))/2 ;
               jcob =abs((nt(2,i+l)-nt(2,i)) /2) ;
               for j=l:quado
                   mark(j) = ycen-gauss(1,j)*jcob;
                   ind = ind+1;
                   gst(1,ind) = nt(1,i);
                   gst(2,ind) = mark(j) ;
                   gst(3,ind) = gauss(2,j) ;
                   gst(4,ind) = jcob;
               end
           end
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  PAGE 236  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%            
%SET UP GAUSS POINTS ALONG DISPLACEMENT BOUNDARY
            gsu = gst
            gsu(1,1:i) = zeros(1,ind);
            qk = zeros(1:2*lthu);
%INTEGRATE FORCES ALONG BOUNDARY
            Imo = (1/12)*Dn3;
            for gt=gst
                gpos = gt(1:2) ;
                weight= gt(3) ;
                jac = gt(4);
                v = domain(gpos,x,dm,numnod) ;
                L = length(v) ;
                en = zeros(1,2*L);
                force=zeros (1,2*L) ;
                [phi,dphix,dphiy] = shape(gpos,dmax,x,v,dm);
                tx = 0 ;
                ty = -(P/(2*Imo))*((D^2)/4-gpos(2,1)^2) ;
                for i=l:L
                     en(2*i-1) = 2*v(i)-1;
                     en(2*i)=2*v(i) ;
                     force(2*i-l)=tx*phi(i);
                     force(2*i) = ty*phi (i) ;
                end
                f(en) = f(en)+jac*weight*force' ;
           end
% INTEGRATE G MATRIX AND Q VECTOR ALONG DISPLACEMENT BOUNDARY
           GG = zeros(numnod*2,lthu*2);
           indl = 0;   ind2=0;
           for i=l:(lthu-l)
               indl=indl+l;
               ml = indl; m2 = ml+l;
               yl = nnu(2,mi) ; y2 = nnu(2,m2) ;
               len = yl-y2;
               for j=l:quado
                   ind2=ind2+1;
                   gpos = gsu(l:2,ind2) ;
                   weight = gsu(3,j);
                   jac = gsu(4,j);
                   facll = (-P*nnu(2,mi)) / (6*young*Imo) ;
                   fac2 = P/ (6*young*Imo) ;
                   xpl = gpos(1,1);
                   ypl = gpos(2,l) ;
                   uxexl = (6*Lb-3*xpl)*xp1+(2+nu)*(yp1^2-(D/2)^2) ;
                   uxexl = uxexl*facll;
                   uyexl = 3*nu*ypln2*(~b-xpl)+0.25*(4+5*nu)*xpl*~n2+(3*Lb-xpl)*xp1n2;
                   uyexl = uyexl*fac2;
 
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  PAGE 237  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%                   
                   
                   Nl = (gpos(2,I)-y2)/len; N2 = 1-N1;
                   qk(2*ml-I) = qk(2*ml-I) -weight* jac*Nl*uxexl ;
                   qk(2*ml) = qk(2*ml) - weight* jac*NI*uyexl ;
                   qk(2*m2-1) = qk(2*m2-1) - weight*jac*N2+uxexl;
                   qk(2*m2) = qk(2*m2) - weight*jac*N2*uyexl;
                   v = domain(gpos ,x ,dm,numnod) ;
                   [phi,dphix,dphiy] = shape(gpos,dmax,x,v,dm);
                   L = length(v);
                   for n=l:L
                       GI = -weight*jac*phi(n)*[Nl 0;0 N1] ;
                       G2 = -weight*jac*phi(n)*[N2 0;0 N2];
                       cl=2*v(n)-l;   c2=2*v(n);   c3=2*m1-l;   c4=2*m1;
                       c5=2*m2-l;c6=2*m2;
                       GG(c1:c2,c3:c4) = GG(c1:c2,c3:c4)+G1;
                       GG(c1:c2,c5:c6) = GG(c1:c2,c5:c6)+ G2 ;
                   end
                end
           end
               
% ENFORCE BC'S USING LAGRANGE MULTIPLIERS
            f = [f;zeros(lthu*2,1)] ;
            f(numnod*2+1:numnod*2+2*lthu,l) = -qk';
            m = sparse([k GG;GG' zeros(lthu*2)]);
            d = m\f ;
            u = d(l:2*numnod) ;
            for i=l:numnod
                u2(1,1) = u(2*i-I) ;
                u2(2,i) = u(2*i);
            end
% SOLVE FOR OUTPUT VARIABLES - DISPLACEMENTS
            for gg=x
                ind = ind+l;
                gpos = gg(1:2);
                v = domain(gpos,x,dm,numod) ;
                [phi,dphix,dphiy] = shape (gpos,dmax,x,v,dm) ;
                displ(2*ind-1) = phi*u2(l,v)' ;
                displ(2*ind) = phi*u2(2,v)' ;
            end
% SOLVE FOR STRESSES AT GAUSS POINTS
            ind = 0;
            for gg=gs
                ind = ind+l;
                gpos = gg(l:2) ;                   
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  PAGE 238  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
 
                weight = gg(3);
                jac = gg(4);
                v = domain(gpos,x,dm,numnod) ;
                L = length(v) ;
                en = zeros (I, 2*L) ;
                [phi, dphix , dphiy] = shape (gpos , dmax , x , v , dm) ;
                Bmat=zeros (3,2*L) ;
                for j=l:L
                    Bmat(l:3, (2*j-I) :2*j) = [dphix(j) 0;O dphiy(j) ;dphiy(j) dphix(j)] ;
                end
                for i=l:L
                    en(2*i-I) = 2*v(i)-I;
                    en(2*i) = 2*v(i);
                end
                stress(l:3,ind) = Dmat*Bmat*u(en);
                stressex(1 ,ind) = (l/Imo)*P*(Lb-gpos(1,I) )*gpos(2,1) ;
                stressex (2, ind) = 0 ;
                stressex(3, ind) = -0.5* (P/Imo) * (0.25*DA2 - gpos (2, I) -2) ;
                err = stress(l:3,ind)-stressex(l:3,ind);
                err2 = weight*jac(0.5*(inv(Dmat)*err)'*(err)) ;
                enorm = enorm + err2;
            end
            enorm = sqrt (enorm)
 
% MESH GENERATION PROGRAM FOR 2D BEAM IN BENDING
       function[x,conn,numcell,numq] = mesh2(length,height,ndivl,ndivw)
% INPUT DATA
            nuncell= ndivw*ndivl;
            numq = (ndivl+l) * (ndivw+l) ;
% SET UP NODAL COORDINATES
           for i = l:(ndivl+l)
               for j = 1 : (ndivw+l)
                   x(l,((ndivw+l)*(i-I) +j))= (length/ndivl)*(i-I);
                   x(2,((ndivw+l)*(i-I) +j))= -(height/ndivw)*(j-l)+height/2;
               end
           end
% SET UP CONNECTIVITY ARRAY
          for j = I :ndivl
              for i = 1:ndivw
                  elemn = (j-l)*ndivw + i;
                  nodet(elemn,l) = elemn + (j-I) ;
                  nodet(elemn,2) = nodet (elemn, I) + 1;
                  nodet(elemn,3) = nodet (elemn, 2) +ndivw+l ;
                  nodet(elemn,4) = nodet (elemn,3) -1 ;
              end
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  PAGE 238  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
 
          end
          conn = nodet';
         
% routine to set up gauss points, jacobian, and weights              
%      function [gs] = egauss(xc,conn,gauss,numcell)
%            index=O ;
%            one = ones(l,4);
%            psiJ = [-l,+l,+l,-1]; etaJ = [-1,-l,+l,+l];
%            l = size(gauss);
%            l = l(2);
%            for e=l:numcell
% % DETERMINE NODES IN EACH CELL
%               for j = 1:4
%                   je=conn(j,e) ;  xe(j)=xc(1,je) ;  ye(j)=xc(2, je) ;
%               end
%               for i=l:l
%                   for j=l:l
%                       index = index+l;
%                       eta=gauss(l ,i);      psi=gauss(l, j) ;
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
%            v(l,l) =-.861136311594052575224;
%            v(1,2) =-.339981043584856264803;
%            v(1,3) = -v(1,2);
%            v(1,4) = -v(l, 1) ;
%            v(2,1) =.347854845137453857373;
%            v(2,2) =.652145154862546142627;
%            v(2,3) = v(2,2);
%            v(2,4) = v(2,I);
% %%%%% ---------------------------------------------------
% %  DETERMINES NODES IN DOMAIN OF INFLUENCE OF POINT GPOS
% %%%%% ---------------------------------------------------
%      function v = domain(gpos,x,dm,nnodes)
%           dif = gpos*ones(l,nnodes)-x;
%           a = dm-abs (dif) ;
%           b = [a;zeros(l ,nnodes)] ;
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
%          t = dm(l:2,v)/dmax;
% % WEIGHTS--W and dw are vectors
%          [w , dwdx , dwdy] = cubwgt (dif , t , v , dmax ,dm) ;
%          B = p.*[w;w;w];
%          pp = zeros(3) ;
%          aa = zeros (3) ;
%          daax = zeros(3);
%          daay = zeros (3) ;
%          for i=l:L
%              pp = p(l:3,i)*p(1:3,i)' ;
%              aa = aa+w(l,i)*pp;
%              daax = daax+dwdx(l,i)*pp;
%              daay = daay+dwdy(l,i)*pp;
%          end
%         pg =[ 1  gpos'] ;
%         [L,U,PERM] = lu(aa);
%         for i=1:3
%             if i==l
%                C = PERM*pg';
%             elseif i==2
%                C = PERM*([0 1 0]' - daax*gam(1:3,1));
%             elseif i==3
%                C = PERM*([0 0 1]' - daay*gam(l:3,1));
%             end
%             Dl = C(1);
%             D2 = (C(2) - L(2,1)*D1) ;
%             D3 = (C(3) - L(3,l)*Dl-L(3,2)*D2) ;
%             gam(3,i) = D3/U(3,3) ;
%             ga1(2,i) = (D2-U(2,3)*gam(3,i))/(U(2,2));
%             gam(1 ,i) =(Dl-U(1,2)*gam(2,i)-U(1,3)*gam(3,i))/U(1,1) ;
%         end
%         phi = gam(l:3,l)'*B;
%         dbx = p.* [ dwdx; dwdx; dwdx] ;
%         dby = p.* [ dwdy; dwdy; dwdy] ;
%         dphix = gam(l:3,2)'*B + gam(l:3,1)'*dbx;
%         dphiy = gam(1:3,3)'*B + gam(l:3,1)'*dby;
% %%%%% ---------------------------------------------------       
%            % CUBIC SPLINE WEIGHT FUNCTION  
% %%%%% --------------------------------------------------- 
% function [w , dwdx , dwdy] = cubwgt (dif , t , v , dmax , dm)
%        1 = length(v) ;
%        for i=l:l
% %%%%    page 241   %%%%%%%%%  
 