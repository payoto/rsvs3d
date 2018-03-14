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
   cout << "volu : index " << index << " | fill " << fill << ", " << 
   target << ", "<< error << " | surfind " << surfind.size();
   for (i=0; unsigned_int(i)<surfind.size();i++){
      cout << "-" << surfind[i];
   }
   cout << endl;
}

void surf::disp() const{
   int i;
   cout << "surf : index " << index << " | fill " << fill << ", " << 
   target << ", "<< error << " | voluind " << voluind.size();
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
   cout << "edge : index " << index <<  " | vertind " << vertind.size();
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
   cout << "vert : index " << index <<  " | edgeind " << edgeind.size();
   for (i=0; unsigned_int(i)<edgeind.size();i++){
      cout << "-" << edgeind[i];
   }
   cout << " | coord " << coord.size();
   for (i=0; unsigned_int(i)<coord.size();i++){
      cout << "-" << coord[i];
   }
   cout << endl;
}

void volu::ChangeIndices(int nVert,int nEdge,int nSurf,int nVolu){
   int i;
   index+=nVolu;
   for (i=0; unsigned_int(i)<surfind.size();i++){
      surfind[i]=surfind[i]+nSurf;
   }
}
void surf::ChangeIndices(int nVert,int nEdge,int nSurf,int nVolu){
   int i;
   index+=nSurf;
   for (i=0; unsigned_int(i)<voluind.size();i++){
      voluind[i]=voluind[i]=nVolu;
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
      surfind[i]=surfind[i]+nSurf;
   }
}
void vert::ChangeIndices(int nVert,int nEdge,int nSurf,int nVolu){
   int i;
   index+=nVert;
   for (i=0; unsigned_int(i)<edgeind.size();i++){
      edgeind[i]=edgeind[i]+nEdge;
   }
}
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
void mesh::PrepareForUse(){
   verts.PrepareForUse();
   edges.PrepareForUse();
   surfs.PrepareForUse();
   volus.PrepareForUse();
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

void mesh::Init(int nVe,int nE, int nS, int nVo){

   verts.Init(nVe);
   edges.Init(nE);
   surfs.Init(nS);
   volus.Init(nVo);

   #ifdef TEST_ARRAYSTRUCTURES
   cout << "Mesh Correctly Assigned!" << endl;
   #endif // TEST_ARRAYSTRUCTURES
}

void mesh::MakeCompatible_inplace(mesh *other) const{
   // Makes other mesh compatible with this to be 
   // merged without index crashes

   int nVert,nEdge,nVolu,nSurf;

   // Define Max indices in current mesh
   this->GetMaxIndex(&nVert,&nEdge,&nVolu,&nSurf);
   other->ChangeIndices(nVert,nEdge,nVolu,nSurf);
}

void mesh::ChangeIndices(int nVert,int nEdge,int nSurf,int nVolu){
   int ii;
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
   }
}

mesh mesh::MakeCompatible(mesh other) const{
   MakeCompatible_inplace(&other);
   return(other);
}

void mesh::Concatenate(const mesh &other){

   this->volus.Concatenate(other.volus);
   this->edges.Concatenate(other.edges);
   this->verts.Concatenate(other.verts);
   this->surfs.Concatenate(other.surfs);
}

void PopulateIndices(mesh *meshin){
   
   meshin->volus.PopulateIndices();
   meshin->edges.PopulateIndices();
   meshin->verts.PopulateIndices();
   meshin->surfs.PopulateIndices();
}

