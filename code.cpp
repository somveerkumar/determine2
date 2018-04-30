#include<stdio.h>
#include<stdlib.h>
#include<math.h>

float determinant(float [][3], float);
void cofactor(float [][3], float);
void transpose(float [][3], float [][3], float);
int n,t;
   double l[25] ,m[25],X[25],Y[25],le[25];
   int ni[25] ,nj[25];
  double K1[8][8], K2[8][8], K3[8][8],k4[8][8] ,global[8][8], C1, C2, C3,c4,A,E,red[3][3],cof[8][8],inv[8][8],Qdisp[8][8], Force[8],L[4];
  float C[4],S[4];
  int main()
  {
   printf("enter the no. of nodes\n");
   scanf("%d",&n);

   printf("enter the no. of elements\n");
   scanf("%d",&t);
   printf("ln element connectivity table\n");
   for(int i=0; i<t; i++)
    {
    scanf("%d\t%d",&ni[i],&nj[i]);  
   }
    printf("element node i node j\n");
    for(int i=0; i<t; i++) 
    {
   printf("%d\t %d\t %d\n",i+1,ni[i],nj[i]);
   }
   printf("\n enter the coordinates of nodes\n");
    for(int i=0;i<n;i++)
   {
    printf("enter the X & Y coordinates of %d\n",i+1);
      scanf("%lf\t%lf",&X[i],&Y[i]);
    }

   printf("Enter the modulus of elasticity in N/mm2\n");
    scanf("%lf",&E);
   printf("Enter the value of area mm2\n");
   scanf("%lf",&A);
    for(int i=0;i<t;i++)
   {
    le[i]=sqrt((pow((X[nj[i]-1]-X[ni[i]-1]),2)+(pow((Y[nj[i]-1]-Y[ni[i]-1]),2))));
   l[i]=(X[nj[i]-1]-X[ni[i]-1])/le[i];
   m[i]=(Y[nj[i]-1]-Y[ni[i]-1])/le[i];
  }
   for(int i=0;i<t;i++)
   {
  
  
    printf("\nle %lf",le[i]);
     printf("\nl%d %lf\t m%d %lf",i,l[i],i,m[i]); 
   printf("\n \n");
 
}
{

   printf("element stiffness matrix 1 \n \n");
}

   float k1[8][8]={{l[0]*l[0],l[0]*m[0],-l[0]*l[0],-l[0]*m[0],0,0,0,0},
  {l[0]*m[0],m[0]*m[0],-l[0]*m[0],-m[0]*m[0],0,0,0,0},
  {-l[0]*l[0],-l[0]*m[0],l[0]*l[0],l[0]*m[0],0,0,0,0},
  {l[0]*m[0],m[0]*m[0],-l[0]*m[0],-m[0]*m[0],0,0,0,0,},
  {0,0,0,0,0,0,0,0}
  ,{0,0,0,0,0,0,0,0},
  {0,0,0,0,0,0,0,0}
  ,{0,0,0,0,0,0,0,0}
  };
  float c1=(A*E)/le[0];
   for(int i=0;i<8;i++)
       {
  for(int j=0;j<8;j++)
  
  	printf(" %f ",c1*k1[i][j]);
  	printf("\n \n");
   }
  
{
   printf("element stiffness matrix 2 \n \n");
}
   float k2[8][8]={{0,0,0,0,0,0,0,0},
   {0,0,0,0,0,0,0,0},
    {0,0,l[1]*l[1],l[1]*m[1],-l[1]*l[1],-l[1]*m[1],0,0}, 
  {0,0,l[1]*m[1],m[1]*m[1],-l[1]*m[1],-m[1]*m[1],0,0},
  {0,0,-l[1]*l[1],-l[1]*m[1],l[1]*l[1],l[1]*m[1],0,0},
  {0,0,l[1]*m[1],m[1]*m[1],-l[1]*m[1],-m[1]*m[1],0,0},{0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0}};

   float c2=(A*E)/le[1];
   
   for(int i=0;i<8;i++)
    {
  for(int j=0;j<8;j++)
  
  	printf(" %f ",c2*k2[i][j]);
  	printf("\n \n");
  
     }
  
  {
   printf("element stiffness matrix 3\n \n");
}
   float k3[8][8]={{0,0,0,0,0,0,0,0},
   {0,0,0,0,0,0,0,0},
    {0,0,l[2]*l[2],l[2]*m[2],0,0,-l[2]*l[2],-l[2]*m[2]}, 
  {0,0,l[2]*m[2],m[2]*m[2],0,0,-l[2]*m[2],-m[2]*m[2]},{0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0},
  {0,0,-l[2]*l[2],-l[2]*m[2],0,0,l[2]*l[2],l[2]*m[2]},
  {0,0,l[2]*m[2],m[2]*m[2],0,0,-l[2]*m[2],-m[2]*m[2]}};

   float c3=(A*E)/le[2];
   
   for(int i=0;i<8;i++)
    {
  for(int j=0;j<8;j++)
  
  	printf(" %f ",c3*k3[i][j]);
  	printf("\n \n");
  
     }
      printf("element stiffness matrix 4 \n \n");
     float k4[8][8]={{0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0},
   {0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0},
    {0,0,0,0,l[3]*l[3],l[3]*m[3],-l[3]*l[3],-l[3]*m[3]}, 
  {0,0,0,0,l[3]*m[3],m[3]*m[3],-l[3]*m[3],-m[3]*m[3]},
  {0,0,0,0,-l[3]*l[3],-l[3]*m[3],l[3]*l[3],l[3]*m[3]},
  {0,0,0,0,l[3]*m[3],m[3]*m[3],-l[3]*m[3],-m[3]*m[3]}};

   float c4=(A*E)/le[3];
   for(int i=0;i<4;i++)
   {
   	L[i]=le[i];
   }
   for(int i=0;i<8;i++)
    {
  for(int j=0;j<8;j++)
  
  	printf(" %f ",c4*k4[i][j]);
  	
  	printf("\n \n");
  
     }
     {
     	printf("Global matrix of 4 element \n \n");
     	
	 }
     
     	  	float sum[8][8];
	
	for(int i=0; i<8; i++)
	{
	
	for(int j=0; j<8; j++)
	{
     sum[i][j]=c1*k1[i][j]+c2*k2[i][j]+c3*k3[i][j]+c4*k4[i][j];	
	}
	}	
	for(int i=0; i<8; i++)
	{
	
	for(int j=0; j<8; j++)
	{
	
	printf("%2f  ",sum[i][j]);
	}	
    printf(" \n  \n ");
	}
	printf("\n reduced global matrix:\n");
	
    float red[3][3]={{sum[2][2],sum[2][4],sum[2][5]},{sum[4][2],sum[4][4],sum[4][5]},{sum[5][2],sum[5][4],sum[5][5]}};
    float a[3][3];
    for(int i=0;i<3;i++)
    {for(int j=0;j<3;j++)
     {
	  printf("%2f \t",red[i][j]);
	  a[i][j]=red[i][j]*.000001;
     }
     printf("\n");
	}
	
	float d;
	 d = determinant(a, 3);
    if (d == 0)
   printf("\nInverse of Entered Matrix is not possible\n");
  else
   cofactor(a,3);
    for(int i=0;i<4;i++)
    {
    	C[i]=l[i];
    	S[i]=m[i];
	}
  }
  void cofactor(float num[3][3], float f)
{
 float b[3][3], fac[3][3];
 int p, q, m, n, i, j;
 for (q = 0;q < f; q++)
 {
   for (p = 0;p < f; p++)
    {
     m = 0;
     n = 0;
     for (i = 0;i < f; i++)
     {
       for (j = 0;j < f; j++)
        {
          if (i != q && j != p)
          {
            b[m][n] = num[i][j];
            if (n < (f - 2))
             n++;
            else
             {
               n = 0;
               m++;
               }
            }
        }
      }
      fac[q][p] = pow(-1, q + p) * determinant(b, f - 1);
    }
  }
  transpose(num, fac, f);
}
float determinant(float a[3][3], float k)
{
  float s = 1, det = 0, b[3][3];
  int i, j, m, n, c;
  if (k == 1)
    {
     return (a[0][0]);
    }
  else
    {
     det = 0;
     for (c = 0; c < k; c++)
       {
        m = 0;
        n = 0;
        for (i = 0;i < k; i++)
          {
            for (j = 0 ;j < k; j++)
              {
                b[i][j] = 0;
                if (i != 0 && j != c)
                 {
                   b[m][n] = a[i][j];
                   if (n < (k - 2))
                    n++;
                   else
                    {
                     n = 0;
                     m++;
                     }
                   }
               }
             }
          det = det + s * (a[0][c] * determinant(b, k - 1));
          s = -1 * s;
          }
    }
 
    return (det);
}
  void transpose(float num[3][3], float fac[3][3], float r)
{
  int i, j;
  float b[3][3], inverse[3][3], d;
 
  for (i = 0;i < r; i++)
    {
     for (j = 0;j < r; j++)
       {
         b[i][j] = fac[j][i];
        }
    }
  d = determinant(num, r);
  for (i = 0;i < r; i++)
    {
     for (j = 0;j < r; j++)
       {
        inverse[i][j] = b[i][j] / d;
        }
    }
   printf("\n\n\nThe inverse of matrix is (*10^-6): \n");
 
   for (i = 0;i < r; i++)
    {
     for (j = 0;j < r; j++)
       {
         printf("\t%f", inverse[i][j]);
         inverse[i][j]=inverse[i][j]*.000001;
        }
    printf("\n");
     }
     float force[3]={0,0,0};
     printf("\n enter the force matrix values:");
     for(i=0;i<3;i++)
     {
	 scanf("%f",&force[i]);
	 printf("\t");
     }
     float dis[3]={0,0,0};
     printf("\n displacement value of node: \n");
     for (i=0;i<3;i++)
     {  
     	dis[i]=inverse[i][0]*force[0]+inverse[i][1]*force[1]+inverse[i][2]*force[2];
     	printf("%f \n",dis[i]);
		 	 }
	 /*float q[8]={0,0,dis[0],0,dis[1],dis[2],0,0};
	 float s1[4]={q[0],q[1],q[2],q[3]};
	 float s2[4]={q[2],q[3],q[4],q[5]};
	 float s3[4]={q[2],q[3],q[6],q[7]};
	 float s4[4]={q[4],q[5],q[6],q[7]};
	 float C1[4]={C[0],S[0],-1*C[0],-1*S[0]};
	 float C2[4]={C[1],S[1],-1*C[1],-1*S[1]};
	 float C3[4]={C[2],S[2],-1*C[2],-1*S[2]};
	 float C4[4]={C[3],S[3],-1*C[3],-1*S[3]};
	 
	 for(i=0;i<4;i++)
	 {
	 	s1[i]=E*s1[i];
	 	s1[i]=s1[i]/le[0];
	 	s2[i]=E*s2[i];
	 	s2[i]=s2[i]/le[1];
	 	s3[i]=E*s3[i];
	 	s3[i]=s3[i]/le[2];
	 	s4[i]=E*s4[i];
	 	s4[i]=s4[i]/le[3];
	 }
	 float stress1,stress2,stress3,stress4;
	 printf("\n stress  values:\n1\n");
	 stress1=C1[0]*s1[0]+C1[1]*s1[1]+C1[2]*s1[2]+C1[3]*s1[3];
	 stress2=C2[0]*s2[0]+C2[1]*s2[1]+C2[2]*s2[2]+C2[3]*s2[3];
	 stress3=C3[0]*s3[0]+C3[1]*s3[1]+C3[2]*s3[2]+C3[3]*s3[3];
	 stress4=C4[0]*s4[0]+C4[1]*s4[1]+C4[2]*s4[2]+C4[3]*s4[3];
	 printf("\n%f",stress1);
	 printf("\n%f",stress2);
	 printf("\n%f",stress3);
	 printf("\n%f",stress4);*/
}
