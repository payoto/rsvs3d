#include <iostream>
#include <stdexcept>
#include <sstream>
#include <functional>

#include "arraystructures.hpp"

// Implementation of features

bool CompareFuncOut(function<void()> func1, function<void()> func2){
   bool compFlag;
   stringstream ss1,ss2;
   auto old_buf = cout.rdbuf(ss1.rdbuf()); 

   func1();
   cout.rdbuf(ss2.rdbuf()); 
   func2();
   std::cout.rdbuf(old_buf);

   compFlag=ss1.str().compare(ss2.str())==0;
   return(compFlag);
   
}


// Class function definitions
// Methods of meshpart : volu surf edge vert
void volu::disp() const{
   int i;
   cout << "volu : index " << index << " | fill " << fill << ", "  <<
   target << ", "<< error << " | isBorder " << isBorder << " | surfind " << surfind.size();
   for (i=0; unsigned_int(i)<surfind.size();i++){
      cout << "-" << surfind[i];
   }
   cout << endl;
}

void surf::disp() const{
   int i;
   cout << "surf : index " << index << " | fill " << fill  << ", " << 
   target << ", "<< error << " | isBorder " << isBorder << " | voluind " << voluind.size();
   for (i=0; unsigned_int(i)<voluind.size();i++){
      cout << "-" << voluind[i];
   }
   cout << " | edgeind " << edgeind.size();
   for (i=0; unsigned_int(i)<edgeind.size();i++){
      cout << "-" << edgeind[i];
   }
   cout << endl;
}

void edge::disp() const{
   int i;
   cout << "edge : index " << index << " | isBorder " << isBorder << " | vertind " << vertind.size();
   for (i=0; unsigned_int(i)<vertind.size();i++){
      cout << "-" << vertind[i];
   }
   cout << " | surfind " << surfind.size();
   for (i=0; unsigned_int(i)<surfind.size();i++){
      cout << "-" << surfind[i];
   }
   cout << endl;
}

void vert::disp() const{
   int i;
   cout << "vert : index " << index << " | isBorder " << isBorder << " | edgeind " << edgeind.size();
   for (i=0; unsigned_int(i)<edgeind.size();i++){
      cout << "-" << edgeind[i];
   }
   cout << " | coord " << coord.size();
   for (i=0; unsigned_int(i)<coord.size();i++){
      cout << "-" << coord[i];
   }
   cout << endl;
}
// Input and output
void volu::write(FILE *fid) const{

   int i;

   fprintf(fid, "%i %.16lf %.16lf %.16lf %i ",index,fill,target, error,int(isBorder));
   fprintf(fid, "%i ",int(surfind.size()));
   for (i=0; unsigned_int(i)<surfind.size();i++){
      fprintf(fid, "%i ",surfind[i]);
   }
   fprintf(fid,"\n");
}

void surf::write(FILE * fid) const{
   int i;

   fprintf(fid, "%i %.16lf %.16lf %.16lf %i ",index,fill,target, error,int(isBorder));
   fprintf(fid, "%i ",int(voluind.size()));
   for (i=0; unsigned_int(i)<voluind.size();i++){
      fprintf(fid, "%i ",voluind[i]);
   }
   fprintf(fid, "%i ",int(edgeind.size()));
   for (i=0; unsigned_int(i)<edgeind.size();i++){
      fprintf(fid, "%i ",edgeind[i]);
   }
   fprintf(fid,"\n");
}

void edge::write(FILE * fid) const{
   int i;

   fprintf(fid, "%i %i ",index,int(isBorder));
   fprintf(fid, "%i ",int(vertind.size()));
   for (i=0; unsigned_int(i)<vertind.size();i++){
      fprintf(fid, "%i ",vertind[i]);
   }
   fprintf(fid, "%i ",int(surfind.size()));
   for (i=0; unsigned_int(i)<surfind.size();i++){
      fprintf(fid, "%i ",surfind[i]);
   }
   fprintf(fid,"\n");
}

void vert::write(FILE * fid) const{
   int i;
   
   fprintf(fid, "%i %i ",index,int(isBorder));
   fprintf(fid, "%i ",int(edgeind.size()));
   for (i=0; unsigned_int(i)<edgeind.size();i++){
      fprintf(fid, "%i ",edgeind[i]);
   }
   fprintf(fid, "%i ",int(coord.size()));
   for (i=0; unsigned_int(i)<coord.size();i++){
      fprintf(fid, "%.16lf ",coord[i]);
   }
   fprintf(fid,"\n");
}

void volu::read(FILE * fid) {

   int i,n;

   fscanf(fid, "%i %lf %lf %lf %i ",&index,&fill,&target, &error, &i);
   isBorder=bool(i);
   fscanf(fid, "%i ",&n);
   surfind.assign(n,0);
   for (i=0; unsigned_int(i)<surfind.size();i++){
      fscanf(fid, "%i ",&surfind[i]);
   }

}

void surf::read(FILE * fid) {
   int i,n;

   fscanf(fid, "%i %lf %lf %lf %i ",&index,&fill,&target, &error, &i);
   isBorder=bool(i);
   fscanf(fid, "%i ",&n);
   voluind.assign(n,0);
   for (i=0; unsigned_int(i)<voluind.size();i++){
      fscanf(fid, "%i ",&voluind[i]);
   }
   fscanf(fid, "%i ",&n);
   edgeind.assign(n,0);
   for (i=0; unsigned_int(i)<edgeind.size();i++){
      fscanf(fid, "%i ",&edgeind[i]);
   }

}

void edge::read(FILE * fid) {
   int i,n;

   fscanf(fid, "%i %i ",&index, &i);
   isBorder=bool(i);
   fscanf(fid, "%i ",&n);
   vertind.assign(n,0);
   for (i=0; unsigned_int(i)<vertind.size();i++){
      fscanf(fid, "%i ",&vertind[i]);
   }
   fscanf(fid, "%i ",&n);
   surfind.assign(n,0);
   for (i=0; unsigned_int(i)<surfind.size();i++){
      fscanf(fid, "%i ",&surfind[i]);
   }

}

void vert::read(FILE * fid) {
   int i,n;
   
   fscanf(fid, "%i %i ",&index, &i);
   isBorder=bool(i);
   fscanf(fid, "%i ",&n);
   edgeind.assign(n,0);
   for (i=0; unsigned_int(i)<edgeind.size();i++){
      fscanf(fid, "%i ",&edgeind[i]);
   }
   fscanf(fid, "%i ",&n);
   coord.assign(n,0);
   for (i=0; unsigned_int(i)<coord.size();i++){
      fscanf(fid, "%lf ",&coord[i]);
   }

}

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-parameter"
void volu::ChangeIndices(int nVert,int nEdge,int nSurf,int nVolu){
   int i;
   index+=nVolu;
   for (i=0; unsigned_int(i)<surfind.size();i++){
      surfind[i]= (surfind[i]>0)? (surfind[i]+nSurf) : surfind[i];
   }
}
void surf::ChangeIndices(int nVert,int nEdge,int nSurf,int nVolu){
   int i;
   index+=nSurf;
   for (i=0; unsigned_int(i)<voluind.size();i++){
      voluind[i]= (voluind[i]>0) ? (voluind[i]+nVolu) : voluind[i];
   }

   for (i=0; unsigned_int(i)<edgeind.size();i++){
      edgeind[i]=edgeind[i]+nEdge;
   }
}
void edge::ChangeIndices(int nVert,int nEdge,int nSurf,int nVolu){
   int i;
   index+=nEdge;
   for (i=0; unsigned_int(i)<vertind.size();i++){
      vertind[i]=vertind[i]+nVert;
   }

   for (i=0; unsigned_int(i)<surfind.size();i++){
      surfind[i]=(surfind[i]>0) ? (surfind[i]+nSurf) : surfind[i];
   }
}
void vert::ChangeIndices(int nVert,int nEdge,int nSurf,int nVolu){
   int i;
   index+=nVert;
   for (i=0; unsigned_int(i)<edgeind.size();i++){
      edgeind[i]=edgeind[i]+nEdge;
   }
}
void mesh::TightenConnectivity(){
   verts.TightenConnectivity();
   surfs.TightenConnectivity();
   edges.TightenConnectivity();
   volus.TightenConnectivity();

}

void mesh::SwitchIndex(int typeInd, int oldInd, int newInd, vector<int> scopeInd)
{

   int ii,jj,kk,newSub,oldSub;
   vector<int> subList;
   HashedVector<int,int> tempSub;
   bool is3DMesh=volus.size()>0;

   if(typeInd==1){
      newSub=verts.find(newInd);
      oldSub=verts.find(oldInd);
      subList=edges.find_list(verts[oldSub].edgeind);
      for (ii=0;ii<int(subList.size());++ii){ // update vertind
         jj=edges(subList[ii])->vertind[1]==oldInd;
         edges[subList[ii]].vertind[jj]=newInd;
         verts[newSub].edgeind.push_back(edges[subList[ii]].index); // update vertex edgeind
         //cout << " " << edges[subList[ii]].index <<  " ";
         for (jj=0;jj<int(verts(oldSub)->edgeind.size());++jj){
            if(verts(oldSub)->edgeind[jj]==edges[subList[ii]].index){
               verts[oldSub].edgeind.erase(
                  verts[oldSub].edgeind.begin()+jj);
               jj--;
            }
         }
      }
      // Hashing has not been invalidated
      edges.isHash=1;
      verts.isHash=1;

   } else if (typeInd==2){
      newSub=edges.find(newInd);

      subList=verts.find_list(edges(edges.find(oldInd))->vertind);
      for (ii=0;ii<int(subList.size());++ii){
         for (jj=0;jj<int(verts(subList[ii])->edgeind.size());++jj){
            if(verts(subList[ii])->edgeind[jj]==oldInd){
               //cout << " " << verts(subList[ii])->index << " " << endl;DisplayVector(verts[subList[ii]].edgeind);cout << endl;
               verts[subList[ii]].edgeind[jj]=newInd;
               //DisplayVector(verts[subList[ii]].edgeind);cout << endl;

            }
         }  
      }
      // Changes the indices of 
      subList=surfs.find_list(edges(edges.find(oldInd))->surfind);
      for (ii=0;ii<int(subList.size());++ii){
         if(subList[ii]!=-1 || is3DMesh){
            for (jj=0;jj<int(surfs(subList[ii])->edgeind.size());++jj){
               if(surfs(subList[ii])->edgeind[jj]==oldInd){
                  surfs[subList[ii]].edgeind[jj]=newInd;
                  edges[newSub].surfind.push_back(surfs[subList[ii]].index);
                  surfs[subList[ii]].isordered=false;
               }
            }
         }  
      }

      edges.isHash=1;
      verts.isHash=1;
      surfs.isHash=1;
      edges.isSetMI=1;
      verts.isSetMI=1;
      surfs.isSetMI=1;


   } else if (typeInd==3){
      newSub=surfs.find(newInd);

      subList=edges.find_list(surfs(surfs.find(oldInd))->edgeind);
      for (ii=0;ii<int(subList.size());++ii){
         for (jj=0;jj<int(edges(subList[ii])->surfind.size());++jj){
            if(edges(subList[ii])->surfind[jj]==oldInd){
               edges[subList[ii]].surfind[jj]=newInd;
               surfs[newSub].edgeind.push_back(edges[subList[ii]].index);
            }
         }  
      }
      surfs.isHash=1;

      subList=volus.find_list(surfs(surfs.find(oldInd))->voluind);
      for (ii=0;ii<int(subList.size());++ii){
         if(subList[ii]!=-1){
            for (jj=0;jj<int(volus(subList[ii])->surfind.size());++jj){
               if(volus(subList[ii])->surfind[jj]==oldInd){
                  volus[subList[ii]].surfind[jj]=newInd;
                  surfs[newSub].voluind.push_back(volus[subList[ii]].index);
               }
            }
         }  
      }

      surfs[newSub].isordered=false;
      surfs.isHash=1;
      edges.isHash=1;
      volus.isHash=1;
      surfs.isSetMI=1;
      edges.isSetMI=1;
      volus.isSetMI=1;

   } else if (typeInd==4){
      newSub=volus.find(newInd);
      subList=surfs.find_list(volus[volus.find(oldInd)].surfind);
      for (ii=0;ii<int(subList.size());++ii){ // update vertind
         //jj=surfs(subList[ii])->voluind[1]==oldInd;
         //surfs[subList[ii]].voluind[jj]=newInd;
         for (jj=0;jj<int(surfs(subList[ii])->voluind.size());++jj){
            if(surfs[subList[ii]].voluind[jj]==oldInd){
               surfs[subList[ii]].voluind[jj]=newInd; 
                  volus[newSub].surfind.push_back(surfs[subList[ii]].index); // update vertex edgeind
               }
            }

         }
      // Hashing has not been invalidated
         volus.isHash=1;
         surfs.isHash=1;
         volus.isSetMI=1;
         surfs.isSetMI=1;
   } else if (typeInd==5){ // Modify vertex index in scoped mode

      newSub=verts.find(newInd);
      oldSub=verts.find(oldInd);

      subList=edges.find_list(scopeInd);
      tempSub.vec=edges.find_list(verts(oldSub)->edgeind);
      tempSub.GenerateHash();
      for (ii=0;ii<int(subList.size());++ii){
         if(tempSub.find(subList[ii])!=-1){
            for (jj=0;jj<int(edges(subList[ii])->vertind.size());++jj){
               if(edges(subList[ii])->vertind[jj]==oldInd){
                  edges[subList[ii]].vertind[jj]=newInd;
                  verts[newSub].edgeind.push_back(edges[subList[ii]].index); 

                  //cout << " " << edges[subList[ii]].index <<  " ";
                  for (kk=0;kk<int(verts(oldSub)->edgeind.size());++kk){
                     if(verts(oldSub)->edgeind[kk]==edges[subList[ii]].index){
                        verts[oldSub].edgeind.erase(
                           verts[oldSub].edgeind.begin()+kk);
                        kk--;
                     }
                  }
               }
            }
         }  
      }
      edges.isHash=1;
      verts.isHash=1;
      edges.isSetMI=1;
      verts.isSetMI=1;

   } else {

      cerr << "Error unknown type of object for index switching" <<endl;
      cerr << "Error in " << __PRETTY_FUNCTION__ << endl;
      throw invalid_argument (" : Type is out of range");
   }

}

void mesh::RemoveIndex(int typeInd, int oldInd)
{

   int ii,jj;
   vector<int> subList;

   if(typeInd==1){

      cerr << "not coded yet" << endl;
      throw;

   } else if (typeInd==2){


      subList=verts.find_list(edges(edges.find(oldInd))->vertind);
      for (ii=0;ii<int(subList.size());++ii){
         for (jj=0;jj<int(verts(subList[ii])->edgeind.size());++jj){
            if(verts(subList[ii])->edgeind[jj]==oldInd){
               verts[subList[ii]].edgeind.erase(
                  verts[subList[ii]].edgeind.begin()+jj);
               jj--;
            }
         }  
      }

      subList=surfs.find_list(edges(edges.find(oldInd))->surfind);
      for (ii=0;ii<int(subList.size());++ii){
         if(subList[ii]!=-1){
            for (jj=0;jj<int(surfs(subList[ii])->edgeind.size());++jj){
               if(surfs(subList[ii])->edgeind[jj]==oldInd){
                  surfs[subList[ii]].edgeind.erase(
                     surfs[subList[ii]].edgeind.begin()+jj);
                  surfs[subList[ii]].isordered=false;
                  jj--;
               }
            }
         }  
      }

      edges.isHash=1;
      verts.isHash=1;
      surfs.isHash=1;
      edges.isSetMI=1;
      verts.isSetMI=1;
      surfs.isSetMI=1;


   } else if (typeInd==3){


      subList=edges.find_list(surfs(surfs.find(oldInd))->edgeind);
      for (ii=0;ii<int(subList.size());++ii){
         for (jj=0;jj<int(edges(subList[ii])->surfind.size());++jj){
            if(edges(subList[ii])->surfind[jj]==oldInd){
               edges[subList[ii]].surfind.erase(
                  edges[subList[ii]].surfind.begin()+jj);
               jj--;
            }
         }  
      }

      subList=volus.find_list(surfs(surfs.find(oldInd))->voluind);
      for (ii=0;ii<int(subList.size());++ii){
         if(subList[ii]!=-1){
            for (jj=0;jj<int(volus(subList[ii])->surfind.size());++jj){
               if(volus(subList[ii])->surfind[jj]==oldInd){
                  volus[subList[ii]].surfind.erase(
                     volus[subList[ii]].surfind.begin()+jj);
                  jj--;
               }
            }
         }  
      }
      
      surfs.isHash=1;
      edges.isHash=1;
      volus.isHash=1;
      surfs.isSetMI=1;
      edges.isSetMI=1;
      volus.isSetMI=1;

   } else if (typeInd==4){

      cerr << "not coded yet" << endl;
      throw;
   } else if (typeInd==5){ // Modify vertex index in scoped mode

      cerr << "not coded yet" << endl;
      throw;
   } else {

      cerr << "Error unknown type of object for index switching" <<endl;
      cerr << "Error in " << __PRETTY_FUNCTION__ << endl;
      throw invalid_argument (" : Type is out of range");
   }

}

void mesh::TestConnectivity(){
   int ii,jj,kk,kk2,errCount, errTot;
   vector<int> testSub;

   errCount=0;
   errTot=0;
   kk=int(verts.size());
   for (ii=0; ii<kk;++ii){
      if(verts(ii)->edgeind.size()==0){
         errCount++;
         cerr << " Test Connectivity Error :" << errCount << " vertex " << verts(ii)->index
               << " Has empty connectivity list "; 
               cerr << endl;

      } else {
         testSub=edges.find_list(verts(ii)->edgeind);
         kk2=testSub.size();
         for(jj=0;jj< kk2; ++jj){
            if (testSub[jj]<0 && verts(ii)->edgeind[jj]!=0){
               errCount++;
               cerr << " Test Connectivity Error :" << errCount << " vertex " << verts(ii)->index
               << " makes unknown reference to edge " << verts(ii)->edgeind[jj] << " list: " ; 
               DisplayVector(verts(ii)->edgeind);
               cerr << endl;
            }
         }
      }
   }
   if (errCount>0){
      cerr << "Test Connectivity vertex (edgeind) Errors :" << errCount << endl;
   }


   errTot+=errCount;
   errCount=0;
   kk=int(edges.size());
   for (ii=0; ii<kk;++ii){
      testSub=verts.find_list(edges(ii)->vertind);
      kk2=testSub.size();
      for(jj=0;jj< kk2; ++jj){
         if (testSub[jj]<0 && edges(ii)->vertind[jj]!=0){
            errCount++;
            cerr << " Test Connectivity Error :" << errCount << " edge " << edges(ii)->index
            << " makes unknown reference to vertex " << edges(ii)->vertind[jj] << " list: " ; DisplayVector(edges(ii)->vertind); 
            cerr << endl;
         }
      }
   }
   if (errCount>0){
      cerr << "Test Connectivity edges (vertind) Errors :" << errCount << endl;
   }


   errTot+=errCount;
   errCount=0;
   for (ii=0; ii<kk;++ii){
      testSub=surfs.find_list(edges(ii)->surfind);
      kk2=testSub.size();
      for(jj=0;jj< kk2; ++jj){
         if (testSub[jj]<0 && edges(ii)->surfind[jj]!=0){
            errCount++;
            cerr << " Test Connectivity Error :" << errCount << " edge " << edges(ii)->index
            << " makes unknown reference to surface " << edges(ii)->surfind[jj] << endl;
         }
      }
   }
   if (errCount>0){
      cerr << "Test Connectivity edges (surfind) Errors :" << errCount << endl;
   }



   errTot+=errCount;
   errCount=0;
   kk=int(surfs.size());
   for (ii=0; ii<kk;++ii){
      testSub=edges.find_list(surfs(ii)->edgeind);
      kk2=testSub.size();
      for(jj=0;jj< kk2; ++jj){
         if (testSub[jj]<0){
            errCount++;
            cerr << " Test Connectivity Error :" << errCount << " surf " << surfs(ii)->index
            << " makes unknown reference to edge " << surfs(ii)->edgeind[jj] << endl;
         }
      }
      if (int(testSub.size())==0){
         errCount++;
         cerr << " Test Connectivity Error :" << errCount << " surf " << surfs(ii)->index
         << " has empty edgeind " <<  endl;
      }
   }
   if (errCount>0) {
      cerr << "Test Connectivity surfs (edgeind) Errors :" << errCount << endl;
   }

   errTot+=errCount;
   errCount=0;
   kk=int(surfs.size());
   for (ii=0; ii<kk;++ii){
      testSub=volus.find_list(surfs(ii)->voluind);
      kk2=testSub.size();
      for(jj=0;jj< kk2; ++jj){
         if (testSub[jj]<0 && surfs(ii)->voluind[jj]!=0){
            errCount++;
            cerr << " Test Connectivity Error :" << errCount << " surf " << surfs(ii)->index
            << " makes unknown reference to volu " << surfs(ii)->voluind[jj] << endl;
         }
      }
   }
   if (errCount>0) {
      cerr << "Test Connectivity surfs (voluind) Errors :" << errCount << endl;
   }


   errTot+=errCount;
   errCount=0;
   kk=int(volus.size());
   for (ii=0; ii<kk;++ii){
      testSub=surfs.find_list(volus(ii)->surfind);
      kk2=testSub.size();
      for(jj=0;jj< kk2; ++jj){
         if (testSub[jj]<0 && volus(ii)->surfind[jj]!=0){
            errCount++;
            cerr << " Test Connectivity Error :" << errCount << " edge " << volus(ii)->index
            << " makes unknown reference to vertex " << volus(ii)->surfind[jj] << endl;
         }
      }
   }
   if (errCount>0){
   cerr << "Test Connectivity volus (surfind) Errors :" << errCount << endl;
   }
   errTot+=errCount;
   if (errTot>0){
      cerr << errTot << "  Total errors were detected in the connectivity list" <<endl;
   }
}

#pragma GCC diagnostic pop
// methods for mesh
void mesh::HashArray(){
   verts.HashArray();
   edges.HashArray();
   surfs.HashArray();
   volus.HashArray();
}
void mesh::SetMaxIndex(){
   verts.SetMaxIndex();
   edges.SetMaxIndex();
   surfs.SetMaxIndex();
   volus.SetMaxIndex();
}
void mesh::SetLastIndex(){
   verts.SetLastIndex();
   edges.SetLastIndex();
   surfs.SetLastIndex();
   volus.SetLastIndex();
}

void mesh::PrepareForUse(bool needOrder){

   verts.isInMesh=true;
   edges.isInMesh=true;
   surfs.isInMesh=true;
   volus.isInMesh=true;


   verts.PrepareForUse();
   edges.PrepareForUse();
   surfs.PrepareForUse();
   volus.PrepareForUse();
   // Additional mesh preparation steps
   if (!borderIsSet){
      this->SetBorders();
   }
   if(needOrder){
      this->OrderEdges();
   }
   verts.ForceArrayReady();
   edges.ForceArrayReady();
   surfs.ForceArrayReady();
   volus.ForceArrayReady();
}

void mesh::GetMaxIndex(int *nVert,int *nEdge,int *nSurf,int *nVolu) const{
   *nVert=verts.GetMaxIndex();
   *nEdge=edges.GetMaxIndex();
   *nSurf=surfs.GetMaxIndex();
   *nVolu=volus.GetMaxIndex();
}
void mesh::disp() const {
   verts.disp();
   edges.disp();
   surfs.disp();
   volus.disp();
}
void mesh::read(FILE *fid) {
   verts.read(fid);
   edges.read(fid);
   surfs.read(fid);
   volus.read(fid);

   verts.isInMesh=true;
   edges.isInMesh=true;
   surfs.isInMesh=true;
   volus.isInMesh=true;
}
void mesh::write(FILE *fid) const {
   //fprintf(fid,"%i \n",int(borderIsSet));
   verts.write(fid);
   edges.write(fid);
   surfs.write(fid);
   volus.write(fid);
}
int mesh::read(const char *str) {
   // Convenience read function taking a single char argument
   FILE *fid;
   fid=fopen(str,"r");
   if (fid!=NULL){
      this->read(fid);
      fclose(fid);
   } else {
      cout << "File " << str << "Could not be opened to read" << endl;
      return(1);
   }
   return(0);
}
int mesh::write(const char *str) const {
   FILE *fid;
   fid=fopen(str,"w");
   if (fid!=NULL){
      this->write(fid);
      fclose(fid);
   } else {
      cout << "File " << str << "Could not be opened to write" << endl;
      return(1);
   }
   return(0);
}
bool mesh::isready() const {
   bool readyforuse=true;
   readyforuse=readyforuse & verts.isready();
   readyforuse=readyforuse & edges.isready();
   readyforuse=readyforuse & surfs.isready();
   readyforuse=readyforuse & volus.isready();

   return(readyforuse);
}



void mesh::displight() const {
   cout << "mesh: vert " << verts.size();
   cout << "; edges " << edges.size();
   cout << "; surfs " << surfs.size();
   cout << "; volus " << volus.size() << endl;
}

void mesh::Init(int nVe,int nE, int nS, int nVo)
{
   borderIsSet=false;

   verts.Init(nVe);
   edges.Init(nE);
   surfs.Init(nS);
   volus.Init(nVo);

   verts.isInMesh=true;
   edges.isInMesh=true;
   surfs.isInMesh=true;
   volus.isInMesh=true;

   #ifdef TEST_ARRAYSTRUCTURES
   cout << "Mesh Correctly Assigned!" << endl;
   #endif // TEST_ARRAYSTRUCTURES
}

void mesh::reserve(int nVe,int nE, int nS, int nVo)
{

   verts.reserve(nVe);
   edges.reserve(nE);
   surfs.reserve(nS);
   volus.reserve(nVo);

   verts.isInMesh=true;
   edges.isInMesh=true;
   surfs.isInMesh=true;
   volus.isInMesh=true;

}

void mesh::MakeCompatible_inplace(mesh &other) const{
   // Makes other mesh compatible with this to be 
   // merged without index crashes

   int nVert,nEdge,nSurf,nVolu;

   // Define Max indices in current mesh
   this->GetMaxIndex(&nVert,&nEdge,&nVolu,&nSurf);
   other.ChangeIndices(nVert,nEdge,nSurf,nVolu);
}

void mesh::ChangeIndices(int nVert,int nEdge,int nSurf,int nVolu){
   /*int ii;
   // Volumes
   for(ii=0;ii<volus.size();++ii){
      volus[ii].ChangeIndices(nVert,nEdge,nSurf,nVolu);
   }
   // Surfaces
   for(ii=0;ii<surfs.size();++ii){
      surfs[ii].ChangeIndices(nVert,nEdge,nSurf,nVolu);
   }
   // Volumes
   for(ii=0;ii<edges.size();++ii){
      edges[ii].ChangeIndices(nVert,nEdge,nSurf,nVolu);
   }
   // Surfaces
   for(ii=0;ii<verts.size();++ii){
      verts[ii].ChangeIndices(nVert,nEdge,nSurf,nVolu);
   }*/
   volus.ChangeIndices(nVert,nEdge,nSurf,nVolu);
   edges.ChangeIndices(nVert,nEdge,nSurf,nVolu);
   surfs.ChangeIndices(nVert,nEdge,nSurf,nVolu);
   verts.ChangeIndices(nVert,nEdge,nSurf,nVolu);
}

mesh mesh::MakeCompatible(mesh other) const{
   MakeCompatible_inplace(other);
   return(other);
}

void mesh::Concatenate(const mesh &other){

   this->volus.Concatenate(other.volus);
   this->edges.Concatenate(other.edges);
   this->verts.Concatenate(other.verts);
   this->surfs.Concatenate(other.surfs);
}

void mesh::PopulateIndices(){

   volus.PopulateIndices();
   edges.PopulateIndices();
   verts.PopulateIndices();
   surfs.PopulateIndices();
}


// Field Specific operations
void surf::OrderEdges(mesh *meshin)
{
   unordered_multimap<int,int> vert2Edge;
   vector<int> edge2Vert,edgeSub,edgeIndOrig;
   int vertCurr,edgeCurr;
   int ii,jj;
   std::pair<unordered_multimap<int,int>::iterator,unordered_multimap<int,int>::iterator> range;
   unordered_multimap<int,int>::iterator it;

   if (edgeind.size()>0){
      sort(edgeind);
      unique(edgeind);

      edgeIndOrig=edgeind;
      edgeSub=meshin->edges.find_list(edgeind);
      edge2Vert=ConcatenateVectorField(meshin->edges,&edge::vertind,edgeSub);

      HashVector(edge2Vert, vert2Edge);

      vertCurr=edge2Vert[0];
      edgeCurr=edgeind[0];
      it=vert2Edge.end();
      for(ii=1;ii<int(edgeind.size());++ii){
         range=vert2Edge.equal_range(vertCurr);
         #ifdef SAFE_ACCESS
         if (range.first==vert2Edge.end()){
            cerr << ii << " vert " << vertCurr << "  ";
            DisplayVector(edge2Vert);
            DisplayVector(edgeind);

            cout << it->second << " " << 1/2 << 2/3 <<  endl;
            throw range_error ("unordered_multimap went beyond its range in OrderEdges");
         }
         if (vert2Edge.count(vertCurr)==1){
            jj=edgeIndOrig[(range.first->second)/2]==edgeCurr;
            DisplayVector(edge2Vert);
            DisplayVector(edgeind);
            cerr <<endl;

            cerr << "Error : Surface does not form closed loop" << endl;
            cerr << "ii is : " << ii << " jj is : " << jj << " count is : " ;
            jj=vert2Edge.count(vertCurr);
            cerr << jj  <<  endl; 
             cerr << "Error in :" << __PRETTY_FUNCTION__ << endl;
            throw range_error ("Found a single vertex - surface is not closed");
         }
         #endif // SAFe_ACCESS
         jj=edgeIndOrig[(range.first->second)/2]==edgeCurr;

         it=range.first;
         if (jj){++it;}

         edgeCurr=edgeIndOrig[(it->second)/2];
         jj=edge2Vert[((it->second)/2)*2]==vertCurr; // Warnign ((x/2)*2) not necessarily equal x
         vertCurr=edge2Vert[((it->second)/2)*2+jj];
         edgeind[ii]=edgeCurr;
      }
      isordered=true;
   }
}

int mesh::OrderEdges(){
   int ii;
   bool kk;

   for (ii = 0; ii < surfs.size(); ++ii)
   {
      if (!surfs(ii)->isready(true)){
         surfs[ii].OrderEdges(this);
         kk=true;
      }
   }
   return(kk);
}

void mesh::SetBorders(){
   int ii,jj,nT;
   if (int(volus.size())>0){
      // Update border status of edges
      for(ii=0;ii<surfs.size();++ii){
         jj=0;
         nT=surfs(ii)->voluind.size();
         surfs[ii].isBorder=false;
         while(jj<nT && !surfs(ii)->isBorder){
            surfs[ii].isBorder=surfs[ii].voluind[jj]==0;
            jj++;
         }
      }
      surfs.ForceArrayReady();

      // Update border status of volus
      for(ii=0;ii<volus.size();++ii){
         jj=0;
         nT=volus(ii)->surfind.size();
         volus[ii].isBorder=false;
         while(jj<nT && !volus(ii)->isBorder){
            volus[ii].isBorder=surfs(surfs.find(volus[ii].surfind[jj]))->isBorder;
            jj++;
         }
      }
      volus.ForceArrayReady();

      // Update border status of edges  
      for(ii=0;ii<edges.size();++ii){
         jj=0;
         nT=edges(ii)->surfind.size();
         edges[ii].isBorder=false;
         while(jj<nT && !edges(ii)->isBorder){
            edges[ii].isBorder=surfs(surfs.find(edges[ii].surfind[jj]))->isBorder;
            jj++;
         }
      }
      edges.ForceArrayReady();

   } else {
      for(ii=0;ii<edges.size();++ii){
         jj=0;
         nT=edges(ii)->surfind.size();
         edges[ii].isBorder=false;
         while(jj<nT && !edges(ii)->isBorder){
            edges[ii].isBorder=(edges[ii].surfind[jj]==0);
            jj++;
         }
      }
      edges.ForceArrayReady();

      for(ii=0;ii<surfs.size();++ii){
         jj=0;
         nT=surfs(ii)->voluind.size();
         surfs[ii].isBorder=false;
         while(jj<nT && !surfs(ii)->isBorder){
            surfs[ii].isBorder=edges(edges.find(surfs[ii].edgeind[jj]))->isBorder;
            jj++;
         }
      }
      surfs.ForceArrayReady();
   }


   // Update border status of edges
   for(ii=0;ii<verts.size();++ii){
      jj=0;
      nT=verts(ii)->edgeind.size();
      verts[ii].isBorder=false;
      while(jj<nT && !verts(ii)->isBorder){
         verts[ii].isBorder=edges(edges.find(verts[ii].edgeind[jj]))->isBorder;
         jj++;
      }
   }
   verts.ForceArrayReady();

   borderIsSet=true;
}

void mesh::ForceCloseContainers(){
   
   int ii,jj,iEdge,iSurf,kk;
   int nVert,nEdge,nSurf,nBlocks;
   bool is3DMesh=volus.size()>0;
   vector<int> vertBlock;


   nBlocks=this->ConnectedVertex(vertBlock);

   nVert=verts.size();
   if (is3DMesh){
      // reassign volumes
      volus.elems.clear();
      volus.Init(nBlocks);
      volus.PopulateIndices();
      volus.HashArray();
      for(ii=0;ii<nVert;ii++){
         nEdge=verts(ii)->edgeind.size();
         for(jj=0;jj<nEdge;++jj){
            iEdge=edges.find(verts(ii)->edgeind[jj]);
            nSurf=edges(iEdge)->surfind.size();
            for (kk=0;kk<nSurf;++kk){
               iSurf=surfs.find(edges(iEdge)->surfind[kk]);
               volus.elems[volus.find(vertBlock[ii])].surfind.push_back(edges(iEdge)->surfind[kk]);
               surfs.elems[iSurf].voluind.clear();
               surfs.elems[iSurf].voluind.push_back(vertBlock[ii]);
               surfs.elems[iSurf].voluind.push_back(0);
            }
         }
      }
   } else {
      // reassign surfaces
      surfs.elems.clear();
      surfs.Init(nBlocks);
      surfs.PopulateIndices();
      surfs.HashArray();

      for(ii=0;ii<nVert;ii++){
         nEdge=verts(ii)->edgeind.size();
         for(jj=0;jj<nEdge;++jj){
            iEdge=edges.find(verts(ii)->edgeind[jj]);
            surfs.elems[surfs.find(vertBlock[ii])].edgeind.push_back(verts(ii)->edgeind[jj]);
            edges.elems[iEdge].surfind.clear();
            edges.elems[iEdge].surfind.push_back(vertBlock[ii]);
            edges.elems[iEdge].surfind.push_back(0);
         }
      }
   }



   verts.ForceArrayReady();
   surfs.ForceArrayReady();
   edges.ForceArrayReady();
   volus.ForceArrayReady();
}

int mesh::ConnectedVertex(vector<int> &vertBlock) const{
   // Fills a vector with a number for each vertex corresponding to a
   // group of connected edges it is part of , can be used close surfaces in 2D or volumes
   // in 3D.
   // Uses a flood fill with queue method


   int nVertExplored,nVerts,nBlocks,nCurr,nEdgesCurr,ii,jj,kk;
   vector<bool> vertStatus; // 1 explored 0 not explored

   vector<int> currQueue, nextQueue; // Current and next queues of indices

   // Preparation of the arrays;
   nVerts=verts.size();
   nBlocks=0;
   nVertExplored=0;

   vertStatus.assign(nVerts,false);
   vertBlock.assign(nVerts,0);
   currQueue.reserve(nVerts/2);
   nextQueue.reserve(nVerts/2);

   
   // While Loop, while not all vertices explored
   while(nVertExplored<nVerts){

      // if currQueue is empty start new block
      if(currQueue.size()<1){

         //cout << "Block " << nBlocks << " - " << nVertExplored << " - " << nVerts << endl;
         ii=0;
         while(vertStatus[ii] && ii<nVerts){
            ii++;
         }
         if (vertStatus[ii]){
            cerr << "Error starting point for loop not found despite max number of vertex not reached" <<endl;
            cerr << "Error in " << __PRETTY_FUNCTION__ << endl;
            throw range_error (" : Starting point for block not found");
         }
         currQueue.push_back(ii);
         nBlocks++;
         
      }
      // Explore current queue
      nCurr=currQueue.size();
      for (ii = 0; ii < nCurr; ++ii){
         if (!vertStatus[currQueue[ii]]){
            vertBlock[currQueue[ii]]=nBlocks;
            nEdgesCurr=verts(currQueue[ii])->edgeind.size();
            for(jj=0;jj<nEdgesCurr;++jj){
               kk=int(edges.isearch(verts(currQueue[ii])->edgeind[jj])->vertind[0]
                  ==verts(currQueue[ii])->index);
               nextQueue.push_back(verts.find(
                  edges.isearch(verts(currQueue[ii])->edgeind[jj])->vertind[kk]));
               #ifdef SAFE_ALGO
               if (verts.find(
                  edges.isearch(verts(currQueue[ii])->edgeind[jj])->vertind[kk])==-1){
                  cerr << "Edge index: " << verts(currQueue[ii])->edgeind[jj] << " vertex index:" <<  
                     edges.isearch(verts(currQueue[ii])->edgeind[jj])->vertind[kk] << endl;
                  cerr << "Edge connected to non existant vertex" <<endl;
                  cerr << "Error in " << __PRETTY_FUNCTION__ << endl;
                  throw range_error (" : Vertex not found");
               }
               
               #endif

            }
            vertStatus[currQueue[ii]]=true;
            nVertExplored++;
         }
      }

      // Reset current queue and set to next queue
      currQueue.clear();
      currQueue.swap(nextQueue);

   }
   return(nBlocks);
}