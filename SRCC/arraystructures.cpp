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
// Input and output
void volu::write(FILE *fid) const{

   int i;

   fprintf(fid, "%i %.16lf %.16lf %.16lf ",index,fill,target, error);
   fprintf(fid, "%i ",int(surfind.size()));
   for (i=0; unsigned_int(i)<surfind.size();i++){
      fprintf(fid, "%i ",surfind[i]);
   }
   fprintf(fid,"\n");
}

void surf::write(FILE * fid) const{
   int i;

   fprintf(fid, "%i %.16lf %.16lf %.16lf ",index,fill,target, error);
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

   fprintf(fid, "%i ",index);
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
   
   fprintf(fid, "%i ",index);
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

   fscanf(fid, "%i %lf %lf %lf ",&index,&fill,&target, &error);
   fscanf(fid, "%i ",&n);
   surfind.assign(n,0);
   for (i=0; unsigned_int(i)<surfind.size();i++){
      fscanf(fid, "%i ",&surfind[i]);
   }

}

void surf::read(FILE * fid) {
   int i,n;

   fscanf(fid, "%i %lf %lf %lf ",&index,&fill,&target, &error);
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

   fscanf(fid, "%i ",&index);
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
   
   fscanf(fid, "%i ",&index);
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
void mesh::read(FILE *fid) {
   verts.read(fid);
   edges.read(fid);
   surfs.read(fid);
   volus.read(fid);
}
void mesh::write(FILE *fid) const {
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

void mesh::Init(int nVe,int nE, int nS, int nVo){

   verts.Init(nVe);
   edges.Init(nE);
   surfs.Init(nS);
   volus.Init(nVo);

   #ifdef TEST_ARRAYSTRUCTURES
   cout << "Mesh Correctly Assigned!" << endl;
   #endif // TEST_ARRAYSTRUCTURES
}

void mesh::MakeCompatible_inplace(mesh &other) const{
   // Makes other mesh compatible with this to be 
   // merged without index crashes

   int nVert,nEdge,nVolu,nSurf;

   // Define Max indices in current mesh
   this->GetMaxIndex(&nVert,&nEdge,&nVolu,&nSurf);
   other.ChangeIndices(nVert,nEdge,nVolu,nSurf);
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
   volus.ChangeIndices(nVert,nEdge,nVolu,nSurf);
   edges.ChangeIndices(nVert,nEdge,nVolu,nSurf);
   surfs.ChangeIndices(nVert,nEdge,nVolu,nSurf);
   verts.ChangeIndices(nVert,nEdge,nVolu,nSurf);
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

