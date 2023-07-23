!Abaqus subroutine for isotropic-isothermal-elasticity with lame parameters in UMAT

SUBROUTINE UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,
      RPL,DDSDDT,DRPLDE,DRPLDT,
      STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,
      NDI,NSHR,NTENS,NSTATV,PROPS,NPROPS,COORDS,DROT,PNEWDT,
      CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,JSTEP,KINC)
C
      INCLUDE 'ABA_PARAM.INC'
C
      CHARACTER*80 CMNAME
      DIMENSION STRESS(NTENS),STATEV(NSTATV),
      DDSDDE(NTENS,NTENS),DDSDDT(NTENS),DRPLDE(NTENS),
      STRAN(NTENS),DSTRAN(NTENS),TIME(2),PREDEF(1),DPRED(1),
      PROPS(NPROPS),COORDS(3),DROT(3,3),DFGRD0(3,3),DFGRD1(3,3),
      JSTEP(4)
!end of base code

!user code for stiffness matrix DDSDDE
  
      real::E,v,G,L1,L2 ; defining variables as decimals datatype
      E=PROPS(1); elasticity
	  v=PROPS(2);poisson's ratio
	  L1=(E*v)/((1.D0+v)*(1-(2.D0*v))); lame parameters
	  L2=E/(2.D0*(1.D0+v)); lame parameters
	  
	   DO i=1, NDI
        DO j=1, NDI
          DDSDDE(i, j)= L1
        END DO
       DDSDDE(i, i)= (2.D0*L2)+L1
       END DO
        
       DO i=NDI+1, NTENS
        DDSDDE(i ,i)= L2
       END DO

!user code to find STRESS matrix
    
	DO i=1, NTENS
     DO j=1, NTENS
      STRESS(j)=STRESS(j) + DDSDDE(j, i) * DSTRAN(i)
     END DO
    END DO
    
	
	RETURN
    END
