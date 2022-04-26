 
 /************************************************************************************************
 *   General swap-based multiple neighborhood tabu search for finding maximum independent set.   *
 *                                                                         	  				     * 
 *   This program is free software; you can redistribute it and/or modify it under the terms     *
 *   of the GNU General Public License as published by the Free Software Foundation; either      *
 *   version 2 of the License, or (at your option) any later version.                            *
 *                                                                         	  				     * 
 *   This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;   *
 *   without even the implied warranty of  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  * 
 *   See the GNU General Public License for more details.                                        *
 * 								                                                                 *
 *   For more information about this algorithm, please visit the website                         *
 *   http://www.info.univ-angers.fr/pub/hao/mis.html or email to: jin@info-univ.angers.fr.       *
 *                                                                                               *
 *                                                                                               *
 ************************************************************************************************/


#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <sstream>
#include <fstream>
#include <string.h>
#include <time.h>
 #include <vector>
#include <ctime>

#define numNb0 1
#define numNb1 2
#define numNb2 3

using namespace std;

char *inFile, *outFile, resultFile[100];
int **edge;
int max_vtx;
int **adjactMatrix, *adjactLen, **adjactIndex;
int *nbAdjactLen;
int *unAdjactLen;
int **beAdjactMatrix, *beAdjactLen, **beAdjactIndex;
int *vertex, *vertexIndex, vertexLen;
int *uVertex, *uVertexIndex, uVertexLen;
int *vertexFlag;
int f, f_best, globe_f;
int *best_vertex, *globe_vertex;
int *tabuListS;
int tabuTenue;

int *nb0, *nb0Index;
int nb0Len;
int *nb1, *nb1Index;
int nb1Len;
int *nb2, *nb2Index;
int nb2Len;

int seed;
int obj;
double dens;
double start_time, end_time, run_time, ave_time, globe_time, ave_best_time;

int eachItersRecord;

void readInitial();
void constructSol();
int multi_tabu_search(int);
void updateAdjVertex(int, int, bool, bool);
void updatesubAdjVertex(int);
void updateInfo(int, int, bool, bool);
void setTabuTenue(int, int, bool, bool);
void clearAll();
bool verify(int, int);

void readInitial()
{
     ifstream FIC;
     FILE *fp;
     FIC.open(inFile);
     if(FIC.fail())
     {
           cout << "### Erreur open, File_Name " << inFile << endl;
           getchar();
           exit(0);                   
     }
     char StrReading[100];
 
     FIC >> StrReading;
     if ( FIC.eof() )
     {
           cout << "### Error open, File_Name " << inFile << endl;
           getchar();
           exit(0);
     }
     int nb_vtx=0, nb_edg=-1, max_edg=0;
     int x1, x2;
     while ( ! FIC.eof() )
     {
           char bidon[50];
           
           if ( strcmp(StrReading, "p" )==0 )
           {
                  FIC >> StrReading;
                  FIC >> max_vtx >> nb_edg;
                  // cout << "Number of vertexes = " << max_vtx << endl;
                  // cout << "Number of edges = " << nb_edg << endl;   
                  
                  edge = new int *[max_vtx];
                  adjactMatrix = new int *[max_vtx];
                  adjactLen = new int [max_vtx];
                  adjactIndex = new int *[max_vtx];
                  nbAdjactLen = new int [max_vtx];
                  unAdjactLen = new int [max_vtx];
                  beAdjactMatrix = new int *[max_vtx];
                  beAdjactLen = new int [max_vtx];
                  beAdjactIndex = new int *[max_vtx];
                  vertex = new int [max_vtx];
                  vertexIndex = new int [max_vtx];
                  vertexFlag = new int [max_vtx];
                  uVertex = new int [max_vtx];
                  uVertexIndex = new int [max_vtx];
                  tabuListS = new int [max_vtx];
                  best_vertex = new int [max_vtx];
                  globe_vertex = new int [max_vtx];
                  
                  nb0 = new int [max_vtx];
                  nb0Index = new int [max_vtx];
                  nb1 = new int [max_vtx];
                  nb1Index = new int [max_vtx];
                  nb2 = new int [max_vtx];
                  nb2Index = new int [max_vtx];                 
                  
                  for(int x=0; x<max_vtx; x++)
                  {
                        edge[x] = new int [max_vtx];
                        adjactMatrix[x] = new int [max_vtx];
                        adjactIndex[x] = new int [max_vtx]; 
                        beAdjactMatrix[x] = new int [max_vtx];  
                        beAdjactIndex[x] = new int [max_vtx];                                             
                  }    
                  
                  for(int x=0; x<max_vtx; x++)
                  {
                          adjactLen[x] = 0;
                          nbAdjactLen[x] = 0;
                          beAdjactLen[x] = 0;
                          vertex[x] = -1;
                          vertexIndex[x] = x;
                          vertexFlag[x] = -1;
                          tabuListS[x] = 0;
                          nb0[x] = x;
                          nb0Index[x] = x;                       
                          
                          for(int y=0; y<max_vtx; y++)
                          {
                                  edge[x][y] = 0;
                                  adjactMatrix[x][y] = 0; 
                                  adjactIndex[x][y] = 0;
                                  beAdjactMatrix[x][y] = 0; 
                                  beAdjactIndex[x][y] = 0;                                   
                          } 
                  } 
                  
           } 

           if ( strcmp(StrReading, "e")==0 )
           {
                       FIC >> x1 >> x2;
                       x1--; x2--;
                       if ( x1<0 || x2<0 || x1>=max_vtx || x2 >=max_vtx )
                       {
                                 cout << "### Error of node : x1="<< x1 << ", x2=" << x2 << endl;
                                 getchar();
                                 exit(0);
                       }
                       edge[x1][x2]=edge[x2][x1]=1;
                       max_edg++;
           }
                 
           FIC >> StrReading;
     }
     
     dens = (float) 2.0 * max_edg/(max_vtx*(max_vtx-1));
     // cout << "Density = " << dens << endl;
     
     if ( 0 && max_edg != nb_edg )
     {
           cout << "### Error de lecture du graphe, nbre aretes : annonce="
                 << nb_edg << ", lu=" << max_edg  << endl;
           getchar();
           exit(0);
     } 
             
     for( int x=0 ; x<max_vtx; x++ )
     {
         adjactLen[x] = 0;
         unAdjactLen[x] = 0;
         for( int y=0; y<max_vtx; y++ )
         {
                 if( edge[x][y] == 1 )
                 {
                    adjactMatrix[x][adjactLen[x]] = y;
                    adjactIndex[x][y] = adjactLen[x];
                    adjactLen[x]++;

                    unAdjactLen[x]++;
                 }
         }
     }
     
     globe_f = 0;
     
     nb0Len = max_vtx;
     nb1Len = 0;
     nb2Len = 0;
	 vertexLen = 0;
     uVertexLen = 0;
     
     eachItersRecord = 0;
                        
     FIC.close();
} 

void constructSol()
{
     int selV;
     while(1)
     {
             if(nb0Len > 0)
             {
                       selV = nb0[rand() % nb0Len];
                       updateInfo(selV, 0, false, false);          
             }
             else
                       break;        
     }     
}

bool expandSetFunc(int Iter)
{
       bool perturb;
       int bestV;
       if(nb0Len != 0)
       {
                 perturb = false;
                 bestV = nb0[rand() % nb0Len];
                 f++;
                 updateInfo(bestV, Iter, perturb, false); 
                 return true;                                  
       }
       else
                 return false;     
}

bool IntensificationStepSpecial(int Iter)
{
       int iterV;
       int best_delt;
       int selV, beV, bestV; 
       bool perturb;
       bool rtnFlag, expandFlag;
       
       best_delt = -1;
       rtnFlag = false;
       
       expandFlag = expandSetFunc(Iter);
       if(expandFlag == true)
               return true;
       else
       {
               for(iterV=0; iterV<nb1Len; iterV++)
               {
                       selV = nb1[iterV];
                       beV = beAdjactMatrix[selV][0];
                           
                       if(nbAdjactLen[beV] > 1)
                       {
                               if(tabuListS[selV] <= Iter)
                               {
                                        if(nbAdjactLen[beV] > best_delt) 
                                        {
                                               best_delt = nbAdjactLen[beV];
                                               bestV = selV; 
                                               rtnFlag = true;                                                     
                                        } 
                                        else if(nbAdjactLen[beV] == best_delt)
                                        {
                                                 if(unAdjactLen[selV] > unAdjactLen[bestV])
                                                 {
                                                        bestV = selV;             
                                                 }
                                                 else if(unAdjactLen[selV] == unAdjactLen[bestV])
                                                 {
                                                         if(rand() % 2)
                                                         {
                                                              bestV = selV;
                                                         }      
                                                 }                               
                                        }       
                               }
                       }                                                            
               }
               if(rtnFlag == true)
               {
                       perturb = false;
                       updateInfo(bestV, Iter, perturb, true);
                       return true;
               }
               else 
                       return false;
       }     
}

bool IntensificationStep(int Iter)
{
       int iterV;
       int best_delt;
       int selV, beV, bestV; 
       bool perturb;
       bool rtnFlag, expandFlag;
       
       best_delt = -1;
       rtnFlag = false;

       expandFlag = expandSetFunc(Iter);
       if(expandFlag == true)
               return true;
       else
       {       
               for(iterV=0; iterV<nb1Len; iterV++)
               {
                       selV = nb1[iterV];
                       beV = beAdjactMatrix[selV][0];
                        
                       if(tabuListS[selV] <= Iter)
                       {
                                if(nbAdjactLen[beV] > best_delt) 
                                {
                                       best_delt = nbAdjactLen[beV];
                                       bestV = selV;  
                                       rtnFlag = true;                                                    
                                } 
                                else if(nbAdjactLen[beV] == best_delt)
                                {
                                         if(unAdjactLen[selV] > unAdjactLen[bestV])
                                         {
                                                bestV = selV;              
                                         }
                                         else if(unAdjactLen[selV] == unAdjactLen[bestV])
                                         {
                                                 if(rand() % 2)
                                                 {
                                                      bestV = selV;
                                                 }      
                                         }                               
                                }     
                       }
               }                                                                          
        
        
               if(rtnFlag == true)
               {
                      perturb = false;
                      updateInfo(bestV, Iter, perturb, false);
                      return true;
               } 
               else
                      return false;
        }        
}

void DiversificationStepSpecial(int Iter)
{
       int bestV, selV, best_delt;
       int x;
       bool perturb;       
       bool selFlag;
       
       perturb = true;
       best_delt = -1;
       selFlag = false;
                                                 
       for(x=0; x<uVertexLen; x++)
       {
               selV = uVertex[x];
               if(unAdjactLen[selV] > best_delt && tabuListS[selV] <= Iter)
               {
                          bestV = selV;
                          best_delt = unAdjactLen[selV];  
                          selFlag = true;        
               }
               else if(unAdjactLen[selV] == best_delt && tabuListS[selV] <= Iter)
               {
                         if(rand() % 2)
                              bestV = selV;
               }
       }
                                                  
       if(selFlag == true)
       {
               f += 1 - beAdjactLen[bestV];                              
               updateInfo(bestV, Iter, perturb, true);
       }        
}

void DiversificationStep(int Iter)
{
       int bestV, selV, best_delt;
       int x;
       bool perturb;       
       bool selFlag;       
       
       perturb = true;
       selFlag = false;
       
       if(rand() % 2)
       {
                   best_delt = -1;
                   for(x=0; x<nb2Len; x++)
                   {
                           selV = nb2[x];
                           if(unAdjactLen[selV] > best_delt && tabuListS[selV] <= Iter)
                           {
                                      bestV = selV;
                                      best_delt = unAdjactLen[selV];
                                      selFlag = true;          
                           }
                           else if(unAdjactLen[selV] == best_delt && tabuListS[selV] <= Iter)
                           {
                                     if(rand() % 2)
                                          bestV = selV;
                           }
                   } 
                   
                   if(selFlag == true)
                   { 
                               f--;                                                                               
                               updateInfo(bestV, Iter, perturb, false);
                   } 
       }
       else
       {
                                                              
                   if(uVertexLen != 0)
                   {
                             bestV = uVertex[rand() % uVertexLen];
                             selFlag = true; 
                   }
                   if(selFlag == true)
                   {
                               f += 1 - beAdjactLen[bestV];                              
                               updateInfo(bestV, Iter, perturb, false);
                   }                          
       }         
}

int multi_tabu_search(int maxIters)
{       
       int Iter;
       bool rtnFlag;
       
       rtnFlag = false;
       
       f = vertexLen;
       f_best = f;
       for(int x=0; x<f; x++)
       {
           best_vertex[x] = vertex[x];
       }
       if(f_best >= obj)
       {
                 return f_best;
       }
       
       for(Iter=1; Iter < maxIters ; Iter++)
       {  
                 if(nb1Len > nb2Len + uVertexLen)
                 {
                           rtnFlag = IntensificationStepSpecial(Iter);
                           if(rtnFlag == false)
                                     DiversificationStepSpecial(Iter);
                 }
                 else
                 {
                           rtnFlag = IntensificationStep(Iter);
                           if(rtnFlag == false)
                                     DiversificationStep(Iter);
                 } 
                 
                 if(f >= f_best)
                 {
                         if(f > f_best)
                         {
                             f_best = f;
                             for(int x=0; x<f_best; x++)
                                    best_vertex[x] = vertex[x];
                         }

                         if(f_best == obj)
                         {
                                 return f_best;
                         }
                 }
            
                 eachItersRecord++;                
       }
       
       return f_best;              
}

void setTabuTenue(int v, int Iter, bool perturbFlag, bool tabuLongFlag)
{
       if(perturbFlag == false)
       {
                  if(tabuLongFlag == false)
                              tabuListS[v] = Iter + tabuTenue + rand() % (nb1Len+1);
                  else
                              tabuListS[v] = Iter + nb1Len; 
       }
       else 
       {
                  tabuListS[v] = Iter + 7; 
       }     
}

void updatesubAdjVertex(int adjV)
{
       int y;
       int subAdjV;
       int beIndex, beV, nbId, nbValue, index, tmpV;
       
       //update the information of beAdjactMatrix
       for(y=0; y<adjactLen[adjV]; y++)
       {
                subAdjV = adjactMatrix[adjV][y];
                
                //update the information of beAdjactMatrix
                if(vertexFlag[subAdjV] != 1)
                {
                        unAdjactLen[adjV]++;
                        unAdjactLen[subAdjV]++;                            
                        
                        beIndex = beAdjactIndex[subAdjV][adjV]; 
                        beAdjactLen[subAdjV]--;

                        beV = beAdjactMatrix[subAdjV][beAdjactLen[subAdjV]];      
                        beAdjactMatrix[subAdjV][beIndex] = beV;
                        beAdjactIndex[subAdjV][beV] = beIndex; 
                        if(beAdjactLen[subAdjV] == 0)
                        {
                                        nb0[nb0Len] = subAdjV;
                                        nb0Index[subAdjV] = nb0Len;
                                        nb0Len++; 
                                        
                                        nb1Len--;
                                        nbId = nb1Index[subAdjV];
                                        nbValue = nb1[nb1Len]; 
                                        nb1[nbId] = nbValue;
                                        nb1Index[nbValue] = nbId;                                                                                                                                                        
                        }
                        else if(beAdjactLen[subAdjV] == numNb0)
                        {
                                        nb1[nb1Len] = subAdjV;
                                        nb1Index[subAdjV] = nb1Len;
                                        nb1Len++; 
                                        
                                        beV = beAdjactMatrix[subAdjV][0];
                                        nbAdjactLen[beV]++;                                            
                                        
                                        nb2Len--;
                                        nbId = nb2Index[subAdjV];
                                        nbValue = nb2[nb2Len]; 
                                        nb2[nbId] = nbValue;
                                        nb2Index[nbValue] = nbId;                                                                                                                                                     
                        } 
                        else if(beAdjactLen[subAdjV] == numNb1)
                        {
                                        nb2[nb2Len] = subAdjV;
                                        nb2Index[subAdjV] = nb2Len;
                                        nb2Len++;
                                        
                                        index = uVertexIndex[subAdjV];
                                        uVertexLen--;
                                        tmpV = uVertex[uVertexLen];
                                        uVertex[index] = tmpV;
                                        uVertexIndex[tmpV] = index;                                                    
                        } 
                                                      
                }                                                     
        }     
}


void updateAdjVertex(int bestV, int Iter, bool perturb, bool tabuLong)
{
       int x;
       int adjV;
       int index, tmpV, nbId, nbValue, beV; 
       bool changeState = false;
       
       for(x=0; x<adjactLen[bestV]; x++)
       {
                adjV = adjactMatrix[bestV][x];

                if(vertexFlag[adjV] != 1)
                {
                        unAdjactLen[adjV]--;
                        changeState = false;
                }
                else if(vertexFlag[adjV] == 1)
                {
                           beAdjactLen[adjV] = 0;
                           
                           index = vertexIndex[adjV];
                           vertexLen--;
                           vertexFlag[adjV] = 0;
                           tmpV = vertex[vertexLen];
                           vertex[index] = tmpV; 
                           vertexIndex[tmpV] = index;
                           changeState = true;                                                   
                           
                           setTabuTenue(adjV, Iter, perturb, tabuLong);
                                      
                           nbAdjactLen[adjV] = 0; 
                           
                           updatesubAdjVertex(adjV);
                }
                     
                //update the information of beAdjactMatrix
                beAdjactMatrix[adjV][beAdjactLen[adjV]] = bestV;
                beAdjactIndex[adjV][bestV] = beAdjactLen[adjV];
                beAdjactLen[adjV]++;
                
                if(beAdjactLen[adjV] == numNb0)
                {
                                 if(changeState == false)
                                 {
                                        nb0Len--;
                                        nbId = nb0Index[adjV];
                                        nbValue = nb0[nb0Len]; 
                                        nb0[nbId] = nbValue;
                                        nb0Index[nbValue] = nbId;
                                 }
                                
                                 nb1[nb1Len] = adjV;
                                 nb1Index[adjV] = nb1Len;
                                 nb1Len++;

                                 nbAdjactLen[bestV]++;
                                                                                                                                         
                }            
                else if(beAdjactLen[adjV] == numNb1)
                {
                                nb1Len--;
                                nbId = nb1Index[adjV];
                                nbValue = nb1[nb1Len]; 
                                nb1[nbId] = nbValue;
                                nb1Index[nbValue] = nbId;
                                
                                beV = beAdjactMatrix[adjV][0];
                                nbAdjactLen[beV]--;
                                
                                nb2[nb2Len] = adjV;
                                nb2Index[adjV] = nb2Len;
                                nb2Len++;                                                                  
                }
                else if(beAdjactLen[adjV] == numNb2)
                {
                                nb2Len--;
                                nbId = nb2Index[adjV];
                                nbValue = nb2[nb2Len]; 
                                nb2[nbId] = nbValue;
                                nb2Index[nbValue] = nbId;

                                uVertex[uVertexLen] = adjV;
                                uVertexIndex[adjV] = uVertexLen;
                                uVertexLen++;                                                                  
                }                        
       
       }     
}

void updateInfo(int bestV, int Iter, bool perturb, bool tabuLong)
{
       int nbId, nbValue;
       
       vertex[vertexLen] = bestV;
       vertexIndex[bestV] = vertexLen;
       vertexLen++;
       vertexFlag[bestV] = 1;
       unAdjactLen[bestV] = 0;
        
       if(beAdjactLen[bestV] == 0)
       {
               nbId = nb0Index[bestV];
               nb0Len--;
               nbValue = nb0[nb0Len]; 
               nb0[nbId] = nbValue;
               nb0Index[nbValue] = nbId;     
       }
       else if(beAdjactLen[bestV] == numNb0)
       {
               nbId = nb1Index[bestV];
               nb1Len--;
               nbValue = nb1[nb1Len]; 
               nb1[nbId] = nbValue;
               nb1Index[nbValue] = nbId;    
       } 
       else if(beAdjactLen[bestV] == numNb1)
       {
               nbId = nb2Index[bestV];
               nb2Len--;
               nbValue = nb2[nb2Len]; 
               nb2[nbId] = nbValue;
               nb2Index[nbValue] = nbId;    
       } 
       else if(beAdjactLen[bestV] >= numNb2)
       {
               nbId = uVertexIndex[bestV];
               uVertexLen--;
               nbValue = uVertex[uVertexLen];
               uVertex[nbId] = nbValue;
               uVertexIndex[nbValue] = nbId;           
       }  
       
       updateAdjVertex(bestV, Iter, perturb, tabuLong);
       
       beAdjactLen[bestV] = 0;                                         
}

void clearAll()
{
          for(int x=0; x<max_vtx; x++)
          {
                  beAdjactLen[x] = 0;
                  nbAdjactLen[x] = 0;
                  unAdjactLen[x] = adjactLen[x];
                  vertex[x] = -1;
                  vertexIndex[x] = x;
                  vertexFlag[x] = -1;
                  nb0[x] = x;
                  nb0Index[x] = x;
                  tabuListS[x] = 0;
                  best_vertex[x] = 0;                  
                                    
                  for(int y=0; y<max_vtx; y++)
                  {
                          beAdjactMatrix[x][y] = 0; 
                          beAdjactIndex[x][y] = 0; 
                  } 
          }  
  
          nb0Len = max_vtx;
          nb1Len = 0;
          nb2Len = 0;
          uVertexLen = 0;
          vertexLen = 0;
}

bool verify(int mSize, int *MIS)
{
      for(int x=0; x<mSize-1; x++)
      {
              for(int y=x+1; y<mSize; y++)
              {
                      if(edge[MIS[x]][MIS[y]] == 1 || MIS[x]==MIS[y])
                           return false;
              }
      }
      return true;
}

void releaseMemory()
{
         
     for(int x=0; x<max_vtx; x++)
	 {
                delete []edge[x];
				delete []adjactMatrix[x];
				delete []adjactIndex[x];
				delete []beAdjactMatrix[x];
				delete []beAdjactIndex[x];
				
	 }
     delete []edge;  
     delete []adjactMatrix;
	 delete []adjactLen;
	 delete []adjactIndex;
	 delete []nbAdjactLen;
	 delete []unAdjactLen;
	 delete []beAdjactMatrix;
	 delete []beAdjactLen;
	 delete []vertex;
	  delete []vertexIndex;
	 delete []vertexFlag;
	 delete []uVertex;
	 delete []uVertexIndex;
	 delete []tabuListS;
	 delete []nb0;
	 delete []nb0Index;
	 delete []nb1;  
     delete []nb1Index;
	 delete []nb2;
	 delete []nb2Index;
	 delete []best_vertex;
	 delete []globe_vertex;
}

int main(int argc, char **argv)
{
       if ( argc == 4 )  
       {
          inFile = argv[1];
          outFile = argv[2];
          obj = atoi(argv[3]);
          tabuTenue = 10;
       }
       else if( argc == 5 )  
       {
          inFile = argv[1];
          outFile = argv[2];
          obj = atoi(argv[3]);
          tabuTenue = atoi(argv[4]);
       }   
       else 
       {
           cout << "Error : the user should input four parameters to run the program." << endl;
           getchar();
           exit(0);
       } 
       
       int maxSize = 0; 
       int runTimes = 5;
       int sucRuns = 0;
       double aveItersRecord = 0.0;
       double aveSize = 0.0;
       int len;
       ave_time = 0.0;
       ave_best_time = 0.0;
       
       len = strlen(outFile);
       strncpy(resultFile,outFile, len-4) ;
       strcat(resultFile, "_r.txt");       

       FILE *fp, *fpRecord;
       bool attainFlag = false;
       
       readInitial();

       int best_size = -1;
       vector<int> best_set;
       
       for(int total_iters=0; total_iters<runTimes; total_iters++)
       {
               seed = (unsigned)time(NULL);
               srand( seed ) ;
               // cout<<"seed ="<<seed<<endl; 
               
               attainFlag = false;
               
               start_time = (double)clock();
               
               eachItersRecord = 0;
               
               for(int iters = 1; iters < 10000; iters++)
               {
               
                   constructSol();
                   
                   maxSize = multi_tabu_search(10000);
                   
                   if(maxSize > globe_f)
                   {
                         globe_f =  maxSize;
                         globe_time = ((double)clock() - start_time) / CLOCKS_PER_SEC;  
                         for(int x=0; x<globe_f; x++)
                                 globe_vertex[x] = best_vertex[x]; // Why is this commented out? For performance reasons? Should we uncomment it if we want to take the best indep.~set?
                         if(globe_f == obj)
                         {
                                    sucRuns++;
                                    attainFlag = true;
                                    break;
                         }
                   }
                   clearAll();
               }
               
               aveSize += globe_f;

               if (globe_f > best_size) {
                best_size = globe_f;
                best_set.resize(best_size);
                for (int i=0; i<globe_f; ++i) 
                  best_set[i] = globe_vertex[i];
               }
               
               end_time = (double)clock();
               run_time = (end_time - start_time) / CLOCKS_PER_SEC;            
               
               if(attainFlag == true)
               {                       
                       ave_time += run_time; 
                       ave_best_time += globe_time;
                       aveItersRecord += eachItersRecord;

                       // fp = fopen(outFile, "a+");
                       // cout<<"maxSize of the independent set is "<<globe_f<<endl;
                       // fprintf(fp, "maxSize of the independent set is %d ", globe_f);
                       // fprintf(fp, " \n ");
                       ///////////
                       // cout<<"best_vertex:";
                       // for(int x=0; x<globe_f; x++)
                       // {
                           // cout<<globe_vertex[x]<<" ";
                           // fprintf(fp, " %d ", globe_vertex[x]);
                       // }
                       // cout<<endl;
                       // fprintf(fp, " \n ");
                       /////////
                       // cout<<"run_time="<<run_time<<endl;
                       // fprintf(fp, "run_time= %f ", run_time);
                       // cout<<"best_run_time="<<globe_time<<endl;
                       // fprintf(fp, "best_run_time= %f ", globe_time);                       
                       // fprintf(fp, " \n ");
                       
                       //if(verify(globe_f, globe_vertex)==true)
                       //{
                       //   cout<<"Correct#####################"<<endl;
                       //   fprintf(fp, "\n Correct#####################\n");
                       //}
                       //else
                       //{
                       //   cout<<"Error&&&&&&&&&&&&&&&&&&&&&&&"<<endl;
                       //   fprintf(fp, "\n Error&&&&&&&&&&&&&&&&&&&&&&&\n");
                       //}
                       
                       // fclose(fp);
                           
                }                   
               else
               {
                       // cout<<"Cannot attain the maximum clique!!!!!!!!!!!!"<<endl; 

                       // fp = fopen(outFile, "a+");
                       // fprintf(fp, "Cannot attain the maximum clique!!!!!!!!!!!!");
                       // cout<<"maxSize of the independent set is "<<globe_f<<endl;
                       // fprintf(fp, "maxSize of the independent set is %d ", globe_f);
                       // fprintf(fp, " \n ");
                       //////////
                       // cout<<"best_vertex:";
                       // for(int x=0; x<globe_f; x++)
                       // {
                           // cout<<globe_vertex[x]<<" ";
                           // fprintf(fp, " %d ", globe_vertex[x]);
                       // }
                       // cout<<endl;
                       // fprintf(fp, " \n ");
                       /////////
                       // cout<<"run_time="<<run_time<<endl;
                       // fprintf(fp, "run_time= %f ", run_time);
                       // cout<<"best_run_time="<<globe_time<<endl;
                       // fprintf(fp, "best_run_time= %f ", globe_time);                        
                       // fprintf(fp, " \n ");
                       
                       //if(verify(globe_f, globe_vertex)==true)
                       //{
                       //   cout<<"Correct#####################"<<endl;
                       //   fprintf(fp, "\n Correct#####################\n");
                       //}
                       //else
                       //{
                       //   cout<<"Error&&&&&&&&&&&&&&&&&&&&&&&"<<endl;
                       //   fprintf(fp, "\n Error&&&&&&&&&&&&&&&&&&&&&&&\n");
                       //}
                     
                       fclose(fp);                                              
               }                    

               clearAll();
               globe_f = 0;
       }
       if(sucRuns != 0)
       {
                  ave_time = ave_time / sucRuns; 
                  ave_best_time = ave_best_time / sucRuns; 
                  aveItersRecord = aveItersRecord / sucRuns;
       }
       else
       {
                  ave_time = 0.0;
                  ave_best_time = 0.0;
                  aveItersRecord = 0.0;
       }
       aveSize = aveSize / runTimes;
       // cout<<"ave_time = " << ave_time <<endl;
       // cout<<"ave_best_time = " << ave_best_time <<endl;
       // fp = fopen(outFile, "a+");
       // fprintf(fp, "\n The success times is %d ", sucRuns);
       // fprintf(fp, "\n The average time is %f ", ave_time);
       // fprintf(fp, "\n The average time of obtaining the best result is %f ", ave_best_time);
       // fprintf(fp, "\n The average iterations is %f ", aveItersRecord);
       // fprintf(fp, "\n The average size of maximum independent set is %f ", aveSize);              
       // fclose(fp);
// 
       // fpRecord = fopen(resultFile, "a+");
       // fprintf(fpRecord, "sucRuns %d", sucRuns);
       // fprintf(fpRecord, " aveTime %f ", ave_time);
       // fprintf(fpRecord, " aveBestTime %f ", ave_best_time);
       // fprintf(fpRecord, " aveIters %f ", aveItersRecord);
       // fprintf(fpRecord, " aveSize %f \n", aveSize);              
       // fclose(fpRecord);     
       for (int i=0; i<best_set.size(); ++i) cout << best_set[i] << " " ;
        cout<<endl;
       
       releaseMemory();	   
       //system("pause"); 
       return 0; 
}
