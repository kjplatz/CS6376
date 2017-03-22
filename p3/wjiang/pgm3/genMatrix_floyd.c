/*
 *generates matrix for floyd programs
 *
 */

#define INFTY ((int) 1 << (sizeof(int) * 8 - 3))

#include<stdio.h>
#include<stdlib.h>

/******************************************************************************/
int main(int argc, char *argv[])
{
  int i, j;
  int n;
  FILE* fp;
  int* Astorage;
  int** A;
  int x;

  if(argc != 3)
    {fprintf(stderr, "\nUsage: genMatrix <n> <outfile>");
      fprintf(stderr, "\nWhere nxn is the size of the matrix into outfile");
      fprintf(stderr, "\n");
      exit(1);
    }

  n = atoi(argv[1]);

  if((fp = fopen(argv[2], "w")) == NULL)
    {fprintf(stderr, "\n*** can't write file %s ***\n", argv[2]);
      exit(1);
    }

  /* write array dimensions n and n */
  fwrite(&n, sizeof(int), 1, fp);
  fwrite(&n, sizeof(int), 1, fp);

  if((Astorage = (int *)malloc(n * n * sizeof(int))) == NULL)
    {fprintf(stderr, "\n*** out of memory ***\n");
      exit(2);
    }

  if((A = (int **)malloc(n * sizeof(int *))) == NULL)
    {fprintf(stderr, "\n*** out of memory ***\n");
      exit(2);
    }

  /* initialize array handle */
  for(i = 0; i < n; ++i)
    A[i] = &Astorage[i * n];

  /* set all values */
  for(i = 0; i < n; ++i)
    for(j = 0; j < n; ++j)
    {
		if(i==j)
			A[i][j]  = 0;
		else if((i-j==1)||(j-i==1)||(i==0&& (((j+1)%(n/8))==0))||(j==0&& (((i+1)%(n/8))==0))  )
			A[i][j]  = 1;
		else
			A[i][j]  = INFTY;
    }

  /* set diagonal to 0 */
  for(i = 0; i < n; ++i)
    A[i][i] = 0;

  /* write to file */
  fwrite(Astorage, sizeof(int), n * n, fp);

  fclose(fp);
  return(0);

}
/******************************************************************************/
/******************************************************************************/
/******************************************************************************/
