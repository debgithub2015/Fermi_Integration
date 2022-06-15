!! This module integrates a quantity with atleast three different method
!! DC 10/20/2014
        PROGRAM Integration

        IMPLICIT NONE

        INTEGER i,j,k,l,n,m,p,q
        INTEGER,dimension(:),allocatable::Natm
        REAL*8,dimension(:),allocatable::tred
        REAL*8 ::one,two,half,pi,ystr,yend,dy,thetaend,thetastr,dtheta,&
                 tmax,tmin,dt,Nmax,Nmin,dNatm,xxint,yyint
        REAL*8,dimension(:,:),allocatable::xxint1,yyint1,&
                                           Refc,Imfc,R2,phi   
       
        ystr=0.d0
        yend=10.d0
        dy=0.025d0
        m=NINT((yend-ystr)/dy)  
           
        one=1.d0
        two=2.d0
        half=0.5d0
        pi=4.d0*atan(1.d0)

        thetastr=0.d0
        thetaend=two*pi
        dtheta=pi*0.0025d0
        n=NINT((thetaend-thetastr)/dtheta)
       
        !WRITE(6,*) n          

        tmin=0.01d0
        tmax=4.d0
        dt=0.025d0
        p=NINT((tmax-tmin)/dt)

        Nmin=2
        Nmax=6
        dNatm=2
        q=(Nmax-Nmin)/dNatm

        allocate(xxint1(1:10000,1:10000),yyint1(1:10000,1:10000),&
                 Refc(1:10000,1:10000),Natm(1:10000),tred(1:10000),&
                 R2(1:10000,1:10000),phi(1:10000,1:10000),&
                 Imfc(1:10000,1:10000))
 
          do k=1,q

            if (k==1) then
               Natm(1)=Nmin
            else
               Natm(k)=Natm(1)+dNatm*(k-1)
            end if

            WRITE(6,'(I5)') Natm(k)

          do l=1,p

            if (l==1) then
               tred(1)=tmin
            else
               tred(l)=tred(1)+dt*(l-1)
            end if
            
          call intsimp2(Natm(k),tred(l),thetastr,thetaend,dtheta,&
                         ystr,yend,dy,xxint,yyint)
          call intsimp22(Natm(k),tred(l),thetastr,thetaend,dtheta,&
                           ystr,yend,dy,xxint,yyint)
          call intsimp32(Natm(k),tred(l),thetastr,thetaend,dtheta,&
                         ystr,yend,dy,xxint,yyint)
          call intsimp382(Natm(k),tred(l),thetastr,thetaend,dtheta,&
                          ystr,yend,dy,xxint,yyint)  
          call intsimp3812(Natm(k),tred(l),thetastr,thetaend,dtheta,&
                           ystr,yend,dy,xxint,yyint)
          call intbodein2(Natm(k),tred(l),thetastr,thetaend,dtheta,&
                           ystr,yend,dy,xxint,yyint)
          
!          write(10,*)"xx",xxint,"yy",yyint           
  
          xxint1(k,l)=(half/pi)*xxint
          yyint1(k,l)=(half/pi)*yyint
           
!          write(10,*) "xxint",xxint1(k,l),"yyint",yyint1(k,l)
 
          IF (ABS(xxint1(k,l)).ge.ABS(yyint1(k,l))) THEN
             R2(k,l)=ABS(xxint1(k,l))*SQRT(1.d0+(yyint1(k,l)/&
                                                xxint1(k,l))**2)
          ELSE IF (ABS(xxint1(k,l)).lt.ABS(yyint1(k,l))) THEN
             R2(k,l)=ABS(yyint1(k,l))*SQRT(1.d0+(xxint1(k,l)/&
                                                yyint1(k,l))**2)
          END IF
           
          phi(k,l)=atan(yyint1(k,l)/xxint1(k,l)) 

          Refc(k,l)=-(1.d0/Natm(k))*Log(R2(k,l))
          Imfc(k,l)=-(phi(k,l)/Natm(k))  
   
          WRITE(6,'(5X,F10.5,X,F10.5,X,F10.5)') tred(l),Refc(k,l),&
                                                           Imfc(k,l)

          end do
          end do 
                    
          deallocate(xxint1,yyint1,Refc,Natm,tred,R2,phi,Imfc)            
 
          END PROGRAM integration

          SUBROUTINE intsimp2(N,t,thetin,thetfin,dthet,yin,yfin,dy,&
                              xxint,yyint)  
               
          IMPLICIT NONE             

          INTEGER, intent(in):: N
          REAL*8, intent(in):: t,thetin,thetfin,dthet,yin,yfin,dy
          REAL*8, intent(out):: xxint,yyint
          REAL*8,external::xx,yy  
          REAL*8 sum1,sum2
          INTEGER::m,i           

          m=NINT((thetfin-thetin)/dthet)

          IF (MOD(m,2).ne.0) m=m+1

          sum1=0.d0            
          sum2=0.d0          

          sum1=sum1+xx(N,t,thetin,yin,yfin,dy)+&
                    xx(N,t,thetfin,yin,yfin,dy)
 
          sum2=sum2+yy(N,t,thetin,yin,yfin,dy)+&
                    yy(N,t,thetfin,yin,yfin,dy)

          do i=1,m-1,2
          sum1=sum1+4.d0*xx(N,t,thetin+i*dthet,yin,yfin,dy)
          sum2=sum2+4.d0*yy(N,t,thetin+i*dthet,yin,yfin,dy)
          end do
 
          do i=2,m-2,2
          sum1=sum1+2.d0*xx(N,t,thetin+i*dthet,yin,yfin,dy)
          sum2=sum2+2.d0*xx(N,t,thetin+i*dthet,yin,yfin,dy)   
          end do

          xxint=(dthet*sum1)/3.d0
          yyint=(dthet*sum2)/3.d0  
 
          !write(10,*) "xxint", xxint,"yyint", yyint 

          END SUBROUTINE       

          SUBROUTINE intsimp22(N,t,thetin,thetfin,dthet,yin,yfin,dy,&
                               xxint,yyint) 

          IMPLICIT NONE             

          INTEGER, intent(in):: N
          REAL*8, intent(in):: t,thetin,thetfin,dthet,yin,yfin,dy
          REAL*8, intent(out):: xxint,yyint
          REAL*8,external::xx,yy  
          REAL*8 sum1,sum2
          INTEGER::m,i           

          m=NINT((thetfin-thetin)/dthet)

          IF (MOD(m,2).ne.0) m=m+1

          sum1=0.d0            
          sum2=0.d0           

          sum1=sum1+xx(N,t,thetin,yin,yfin,dy)+&
                    xx(N,t,thetfin,yin,yfin,dy)
 
          sum2=sum2+yy(N,t,thetin,yin,yfin,dy)+&
                    yy(N,t,thetfin,yin,yfin,dy)

          do i=1,2*m-1,2
          sum1=sum1+4.d0*xx(N,t,thetin+i*dthet/2.d0,yin,yfin,dy)
          sum2=sum2+4.d0*yy(N,t,thetin+i*dthet/2.d0,yin,yfin,dy)
          end do
 
          do i=2,2*m-2,2
          sum1=sum1+2.d0*xx(N,t,thetin+i*dthet/2.d0,yin,yfin,dy)
          sum2=sum2+2.d0*xx(N,t,thetin+i*dthet/2.d0,yin,yfin,dy)   
          end do

          xxint=(dthet*sum1)/3.d0/2.d0
          yyint=(dthet*sum2)/3.d0/2.d0

!          write(10,*) "xxint", xxint,"yyint", yyint

          END SUBROUTINE  
           
          SUBROUTINE intsimp32(N,t,thetin,thetfin,dthet,yin,yfin,dy,&
                               xxint,yyint)     

          IMPLICIT NONE         

          INTEGER, intent(in):: N
          REAL*8, intent(in):: t,thetin,thetfin,dthet,yin,yfin,dy
          REAL*8, intent(out):: xxint,yyint
          REAL*8,external::xx,yy  
          REAL*8 sum1,sum2
          INTEGER::m,i           

          m=NINT((thetfin-thetin)/dthet)

          IF (MOD(m,2).ne.0) m=m+1

          sum1=0.d0            
          sum2=0.d0            

          sum1=sum1+(3.d0/8.d0)*(xx(N,t,thetin,yin,yfin,dy)+&
                               xx(N,t,thetfin,yin,yfin,dy))
          sum1=sum1+(7.d0/6.d0)*(xx(N,t,thetin+dthet,yin,yfin,dy)+&
                               xx(N,t,thetfin-dthet,yin,yfin,dy)) 
          sum1=sum1+(23.d0/24.d0)*(xx(N,t,thetin+2.d0*dthet,yin,yfin,dy)+&
                               xx(N,t,thetfin-2.d0*dthet,yin,yfin,dy))
 
          sum2=sum2+3.d0/8.d0*(yy(N,t,thetin,yin,yfin,dy)+&
                               yy(N,t,thetfin,yin,yfin,dy))
          sum2=sum2+(7.d0/6.d0)*(yy(N,t,thetin+dthet,yin,yfin,dy)+&
                               yy(N,t,thetfin-dthet,yin,yfin,dy))
          sum2=sum2+(23.d0/24.d0)*(yy(N,t,thetin+2.d0*dthet,yin,yfin,dy)+&
                               yy(N,t,thetfin-2.d0*dthet,yin,yfin,dy))

          do i=3,m-3
          sum1=sum1+xx(N,t,thetin+i*dthet,yin,yfin,dy)
          sum2=sum2+yy(N,t,thetin+i*dthet,yin,yfin,dy)
          end do
 
          xxint=(dthet*sum1)
          yyint=(dthet*sum2)      

!          write(10,*) "xxint", xxint,"yyint", yyint         
 
          END SUBROUTINE 

          SUBROUTINE intsimp382(N,t,thetin,thetfin,dthet,yin,yfin,dy,&
                                xxint,yyint)
          
          IMPLICIT NONE
   
          INTEGER, intent(in):: N
          REAL*8, intent(in):: t,thetin,thetfin,dthet,yin,yfin,dy
          REAL*8, intent(out):: xxint,yyint
          REAL*8,external::xx,yy  
          REAL*8 sum1,sum2
          INTEGER::m,i           

          m=NINT((thetfin-thetin)/dthet)

          IF (MOD(m,2).ne.1) m=m+1

          sum1=0.d0            
          sum2=0.d0  

          sum1=sum1+(xx(N,t,thetin,yin,yfin,dy)+ &
                     xx(N,t,thetfin,yin,yfin,dy))
          sum2=sum2+(yy(N,t,thetin,yin,yfin,dy)+ &
                     yy(N,t,thetfin,yin,yfin,dy))

          do i=1,m-1,3
          sum1=sum1+3.d0*xx(N,t,thetin+i*dthet,yin,yfin,dy)
          sum2=sum2+3.d0*yy(N,t,thetin+i*dthet,yin,yfin,dy)
          end do
          
          do i=2,m-2,3
          sum1=sum1+3.d0*xx(N,t,thetin+i*dthet,yin,yfin,dy)
          sum2=sum2+3.d0*yy(N,t,thetin+i*dthet,yin,yfin,dy)
          end do
 
          do i=3,m-3,3
          sum1=sum1+2.d0*xx(N,t,thetin+i*dthet,yin,yfin,dy)
          sum2=sum2+2.d0*yy(N,t,thetin+i*dthet,yin,yfin,dy)  
          end do

          xxint=sum1*dthet*3.d0/8.d0 
          yyint=sum2*dthet*3.d0/8.d0 

!          write(10,*) "xxint", xxint,"yyint", yyint

          END SUBROUTINE
         
          SUBROUTINE intsimp3812(N,t,thetin,thetfin,dthet,yin,yfin,dy,&
                                xxint,yyint)
          
          IMPLICIT NONE
   
          INTEGER, intent(in):: N
          REAL*8, intent(in):: t,thetin,thetfin,dthet,yin,yfin,dy
          REAL*8, intent(out):: xxint,yyint
          REAL*8,external::xx,yy  
          REAL*8 sum1,sum2
          INTEGER::m,i           

          m=NINT((thetfin-thetin)/dthet)

          IF (MOD(m,2).ne.1) m=m+1

          sum1=0.d0            
          sum2=0.d0  

          sum1=sum1+(xx(N,t,thetin,yin,yfin,dy)+ &
                     xx(N,t,thetfin,yin,yfin,dy))
          sum2=sum2+(yy(N,t,thetin,yin,yfin,dy)+ &
                     yy(N,t,thetfin,yin,yfin,dy))

          do i=1,4*m-1,3
          sum1=sum1+3.d0*xx(N,t,thetin+i*dthet/4.d0,yin,yfin,dy)
          sum2=sum2+3.d0*yy(N,t,thetin+i*dthet/4.d0,yin,yfin,dy)
          end do
          
          do i=2,4*m-2,3
          sum1=sum1+3.d0*xx(N,t,thetin+i*dthet/4.d0,yin,yfin,dy)
          sum2=sum2+3.d0*yy(N,t,thetin+i*dthet/4.d0,yin,yfin,dy)
          end do
 
          do i=3,4*m-3,3
          sum1=sum1+2.d0*xx(N,t,thetin+i*dthet/4.d0,yin,yfin,dy)
          sum2=sum2+2.d0*yy(N,t,thetin+i*dthet/4.d0,yin,yfin,dy)  
          end do

          xxint=sum1*dthet*3.d0/8.d0/4.d0 
          yyint=sum2*dthet*3.d0/8.d0/4.d0 

!          write(10,*) "xxint", xxint,"yyint", yyint

          END SUBROUTINE

          SUBROUTINE intbodein2(N,t,thetin,thetfin,dthet,yin,yfin,dy,&
                                xxint,yyint)
          
          IMPLICIT NONE
   
          INTEGER, intent(in):: N
          REAL*8, intent(in):: t,thetin,thetfin,dthet,yin,yfin,dy
          REAL*8, intent(out):: xxint,yyint
          REAL*8,external::xx,yy  
          REAL*8 sum1,sum2
          INTEGER::m,i           

          m=NINT((thetfin-thetin)/dthet)

          IF (MOD(m,2).ne.0) m=m+1

          sum1=0.d0            
          sum2=0.d0  

          sum1=sum1+14.d0*(xx(N,t,thetin,yin,yfin,dy)+ &
                     xx(N,t,thetfin,yin,yfin,dy))
          sum2=sum2+14.d0*(yy(N,t,thetin,yin,yfin,dy)+ &
                     yy(N,t,thetfin,yin,yfin,dy))

          do i=1,5*m-1,4
          sum1=sum1+64.d0*xx(N,t,thetin+i*dthet/5.d0,yin,yfin,dy)
          sum2=sum2+64.d0*yy(N,t,thetin+i*dthet/5.d0,yin,yfin,dy)
          end do
          
          do i=2,5*m-2,4
          sum1=sum1+24.d0*xx(N,t,thetin+i*dthet/5.d0,yin,yfin,dy)
          sum2=sum2+24.d0*yy(N,t,thetin+i*dthet/5.d0,yin,yfin,dy)
          end do
 
          do i=3,5*m-3,4
          sum1=sum1+64.d0*xx(N,t,thetin+i*dthet/5.d0,yin,yfin,dy)
          sum2=sum2+64.d0*yy(N,t,thetin+i*dthet/5.d0,yin,yfin,dy)  
          end do

          do i=4,5*m-3,4
          sum1=sum1+2.d0*14.d0*xx(N,t,thetin+i*dthet/5.d0,yin,yfin,dy)
          sum2=sum2+2.d0*14.d0*yy(N,t,thetin+i*dthet/5.d0,yin,yfin,dy)  
          end do          

          xxint=sum1*dthet/45.d0/5.d0 
          yyint=sum2*dthet/45.d0/5.d0 

!          write(10,*) "xxint", xxint,"yyint", yyint

          END SUBROUTINE

          REAL*8 FUNCTION xx(N,t,theta,yin,yfin,dy) 

          IMPLICIT NONE
          INTEGER,intent(in)::N
          REAL*8,intent(in)::t,theta,yin,yfin,dy
          REAL*8 t32,theta1,Nt32gR,Nt32gI,gRint,gIint

          t32=t**(3.d0/2.d0)
          
         call intsimp1(yin,yfin,dy,theta,gRint,gIint) 
         call intsimp21(yin,yfin,dy,theta,gRint,gIint) 
         call intsimp31(yin,yfin,dy,theta,gRint,gIint)
         call intsimp381(yin,yfin,dy,theta,gRint,gIint)
         call intsimp3811(yin,yfin,dy,theta,gRint,gIint)                 
         call intbodein1(yin,yfin,dy,theta,gRint,gIint)  
         
          Nt32gR=N*t32*gRint
          Nt32gI=N*t32*gIint
          theta1=N*theta+Nt32gI
           
          xx=exp(Nt32gR)*cos(theta1)
        
          write(12,*) "theta",theta,"gRint",gRint,"gIint",gIint
!          write(10,*) "theta",theta,"xx", xx
            
          END  FUNCTION  
            
          REAL*8 FUNCTION yy(N,t,theta,yin,yfin,dy) 

          IMPLICIT NONE
          INTEGER,intent(in)::N
          REAL*8,intent(in)::t,theta,yin,yfin,dy
          REAL*8 t32,theta1,Nt32gR,Nt32gI,gRint,gIint

          t32=t**(3.d0/2.d0)
          
         call intsimp1(yin,yfin,dy,theta,gRint,gIint) 
         call intsimp21(yin,yfin,dy,theta,gRint,gIint) 
         call intsimp31(yin,yfin,dy,theta,gRint,gIint)
         call intsimp381(yin,yfin,dy,theta,gRint,gIint)
         call intsimp3811(yin,yfin,dy,theta,gRint,gIint)                 
         call intbodein1(yin,yfin,dy,theta,gRint,gIint) 
 
          Nt32gR=N*t32*gRint
          Nt32gI=N*t32*gIint
          theta1=N*theta+Nt32gI
           
          yy=exp(Nt32gR)*sin(theta1)
           
!          write(10,*) "theta",theta,"yy",yy
            
          END FUNCTION  

          SUBROUTINE intsimp1(yin,yfin,dy,theta,gRint,gIint)
           
          IMPLICIT NONE
          REAL*8,intent(in)::yin,yfin,dy,theta
          REAL*8,intent(out)::gRint,gIint
          REAL*8,external::gR,gI  
          REAL*8 sum1,sum2
          INTEGER::m,i  
           
          m=NINT((yfin-yin)/dy)

          IF (MOD(m,2).ne.0) m=m+1

          sum1=0.d0            
          sum2=0.d0          

          sum1=sum1+gR(theta,yin)+gR(theta,yfin)
 
          sum2=sum2+gI(theta,yin)+gI(theta,yfin)

          do i=1,m-1,2
          sum1=sum1+4.d0*gR(theta,yin+i*dy)
          sum2=sum2+4.d0*gI(theta,yin+i*dy)
          end do
 
          do i=2,m-2,2
          sum1=sum1+2.d0*gR(theta,yin+i*dy)
          sum2=sum2+2.d0*gI(theta,yin+i*dy)   
          end do

          gRint=(dy*sum1)/3.d0
          gIint=(dy*sum2)/3.d0

          !WRITE(10,*) "theta",theta
          !WRITE(10,*) "gRint",gRint
          !WRITE(10,*) "gIint",gIint  
                                 
          END SUBROUTINE
           
          SUBROUTINE intsimp21(yin,yfin,dy,theta,gRint,gIint)
           
          IMPLICIT NONE
          REAL*8,intent(in)::yin,yfin,dy,theta
          REAL*8,intent(out)::gRint,gIint
          REAL*8,external::gR,gI  
          REAL*8 sum1,sum2
          INTEGER::m,i  
           
          m=NINT((yfin-yin)/dy)

          IF (MOD(m,2).ne.0) m=m+1

          sum1=0.d0            
          sum2=0.d0          

          sum1=sum1+gR(theta,yin)+gR(theta,yfin)
 
          sum2=sum2+gI(theta,yin)+gI(theta,yfin)

          do i=1,2*m-1,2
          sum1=sum1+4.d0*gR(theta,yin+i*dy/2.d0)
          sum2=sum2+4.d0*gI(theta,yin+i*dy/2.d0)
          end do
 
          do i=2,2*m-2,2
          sum1=sum1+2.d0*gR(theta,yin+i*dy/2.d0)
          sum2=sum2+2.d0*gI(theta,yin+i*dy/2.d0)   
          end do

          gRint=(dy*sum1)/3.d0/2.d0
          gIint=(dy*sum2)/3.d0/2.d0
                                 
!          WRITE(10,*) "theta",theta
!          WRITE(10,*) "gRint",gRint
!          WRITE(10,*) "gIint",gIint          

          END SUBROUTINE

          SUBROUTINE intsimp31(yin,yfin,dy,theta,gRint,gIint)
           
          IMPLICIT NONE
          REAL*8,intent(in)::yin,yfin,dy,theta
          REAL*8,intent(out)::gRint,gIint
          REAL*8,external::gR,gI  
          REAL*8 sum1,sum2
          INTEGER::m,i  
           
          m=NINT((yfin-yin)/dy)

          IF (MOD(m,2).ne.0) m=m+1

          sum1=0.d0            
          sum2=0.d0          

          sum1=sum1+(3.d0/8.d0)*(gR(theta,yin)+gR(theta,yfin))
          sum1=sum1+(7.d0/6.d0)*(gR(theta,yin+dy)+gR(theta,yfin-dy))
          sum1=sum1+(23.d0/24.d0)*(gR(theta,yin+2.d0*dy)+gR(theta,yfin-2.d0*dy)) 
 
          sum2=sum2+(3.d0/8.d0)*(gI(theta,yin)+gI(theta,yfin))
          sum2=sum2+(7.d0/6.d0)*(gI(theta,yin+dy)+gI(theta,yfin-dy))
          sum2=sum2+(23.d0/24.d0)*(gI(theta,yin+2.d0*dy)+gI(theta,yfin-2.d0*dy))

          do i=3,m-3
          sum1=sum1+gR(theta,yin+i*dy)
          sum2=sum2+gI(theta,yin+i*dy)
          end do

          gRint=(dy*sum1)
          gIint=(dy*sum2)
          
!         WRITE(10,*) "theta",theta
!         WRITE(10,*) "gRint",gRint
!         WRITE(10,*) "gIint",gIint        
                       
          END SUBROUTINE
            
          SUBROUTINE intsimp381(yin,yfin,dy,theta,gRint,gIint)
            
          IMPLICIT NONE
          REAL*8,intent(in)::yin,yfin,dy,theta
          REAL*8,intent(out)::gRint,gIint
          REAL*8,external::gR,gI  
          REAL*8 sum1,sum2
          INTEGER::m,i  
           
          m=NINT((yfin-yin)/dy)

          IF (MOD(m,2).ne.1) m=m+1

          sum1=0.d0            
          sum2=0.d0

          sum1=sum1+(gR(theta,yin)+gR(theta,yfin))
          sum2=sum2+(gI(theta,yin)+gI(theta,yfin))

          do i=1,m-1,3
          sum1=sum1+3.d0*gR(theta,yin+i*dy)
          sum2=sum2+3.d0*gI(theta,yin+i*dy)
          end do
          
          do i=2,m-2,3
          sum1=sum1+3.d0*gR(theta,yin+i*dy)
          sum2=sum2+3.d0*gI(theta,yin+i*dy)
          end do
 
          do i=3,m-3,3
          sum1=sum1+2.d0*gR(theta,yin+i*dy)
          sum2=sum2+2.d0*gI(theta,yin+i*dy)  
          end do

          gRint=sum1*dy*3.d0/8.d0 
          gIint=sum2*dy*3.d0/8.d0        
  
!          WRITE(10,*) "theta",theta
!          WRITE(10,*) "gRint",gRint
!          WRITE(10,*) "gIint",gIint

          END SUBROUTINE

          SUBROUTINE intsimp3811(yin,yfin,dy,theta,gRint,gIint)
            
          IMPLICIT NONE
          REAL*8,intent(in)::yin,yfin,dy,theta
          REAL*8,intent(out)::gRint,gIint
          REAL*8,external::gR,gI  
          REAL*8 sum1,sum2
          INTEGER::m,i  
           
          m=NINT((yfin-yin)/dy)

          IF (MOD(m,2).ne.1) m=m+1

          sum1=0.d0            
          sum2=0.d0

          sum1=sum1+(gR(theta,yin)+gR(theta,yfin))
          sum2=sum2+(gI(theta,yin)+gI(theta,yfin))

          do i=1,4*m-1,3
          sum1=sum1+3.d0*gR(theta,yin+i*dy/4.d0)
          sum2=sum2+3.d0*gI(theta,yin+i*dy/4.d0)
          end do
          
          do i=2,4*m-2,3
          sum1=sum1+3.d0*gR(theta,yin+i*dy/4.d0)
          sum2=sum2+3.d0*gI(theta,yin+i*dy/4.d0)
          end do
 
          do i=3,4*m-3,3
          sum1=sum1+2.d0*gR(theta,yin+i*dy/4.d0)
          sum2=sum2+2.d0*gI(theta,yin+i*dy/4.d0)  
          end do

          gRint=sum1*dy*3.d0/8.d0/4.d0 
          gIint=sum2*dy*3.d0/8.d0/4.d0        
  
!         WRITE(10,*) "theta",theta
!         WRITE(10,*) "gRint",gRint
!         WRITE(10,*) "gIint",gIint

          END SUBROUTINE

          SUBROUTINE intbodein1(yin,yfin,dy,theta,gRint,gIint)
            
          IMPLICIT NONE
          REAL*8,intent(in)::yin,yfin,dy,theta
          REAL*8,intent(out)::gRint,gIint
          REAL*8,external::gR,gI  
          REAL*8 sum1,sum2
          INTEGER::m,i  
           
          m=NINT((yfin-yin)/dy)

          IF (MOD(m,2).ne.0) m=m+1

          sum1=0.d0            
          sum2=0.d0

          sum1=sum1+14.d0*(gR(theta,yin)+gR(theta,yfin))
          sum2=sum2+14.d0*(gI(theta,yin)+gI(theta,yfin))

          do i=1,5*m-1,4
          sum1=sum1+64.d0*gR(theta,yin+i*dy/5.d0)
          sum2=sum2+64.d0*gI(theta,yin+i*dy/5.d0)
          end do
          
          do i=2,5*m-2,4
          sum1=sum1+24.d0*gR(theta,yin+i*dy/5.d0)
          sum2=sum2+24.d0*gI(theta,yin+i*dy/5.d0)
          end do
 
          do i=3,5*m-3,4
          sum1=sum1+64.d0*gR(theta,yin+i*dy/5.d0)
          sum2=sum2+64.d0*gI(theta,yin+i*dy/5.d0)  
          end do

          do i=4,5*m-4,4
          sum1=sum1+2.d0*14.d0*gR(theta,yin+i*dy/5.d0)
          sum2=sum2+2.d0*14.d0*gI(theta,yin+i*dy/5.d0)  
          end do
          
          gRint=sum1*dy/45.d0/5.d0
          gIint=sum2*dy/45.d0/5.d0        
  
!          WRITE(10,*) "theta",theta
!          WRITE(10,*) "gRint",gRint
!          WRITE(10,*) "gIint",gIint

          END SUBROUTINE
           
          REAL*8 FUNCTION gR(theta,y)
           
          IMPLICIT NONE
          REAL*8,intent(in)::theta,y
          REAL*8 :: gRnom,gRdenom
          REAL*8 :: y32,expy,expysin,expycos
           
          y32=y**(3.d0/2.d0)
          expy=exp(y)
          expysin=expy*sin(theta)
          expycos=expy*cos(theta)

          gRnom=y32*(expycos+1.d0)
          gRdenom=(expycos+1.d0)**2+(expysin)**2            
     
          gR=gRnom/gRdenom

!          WRITE(10,*) "theta",theta,"gR", gR

          END FUNCTION 

          REAL*8 FUNCTION gI(theta,y)
           
          IMPLICIT NONE
          REAL*8,intent(in)::theta,y
          REAL*8 :: gInom,gIdenom
          REAL*8 :: y32,expy,expysin,expycos
           
          y32=y**(3.d0/2.d0)
          expy=exp(y)
          expysin=expy*sin(theta)
          expycos=expy*cos(theta)

          gInom=y32*expysin
          gIdenom=(expycos+1.d0)**2+(expysin)**2            
     
          gI=-gInom/gIdenom

!          WRITE(10,*) "theta",theta,"gI",gI  
    
          END FUNCTION 
