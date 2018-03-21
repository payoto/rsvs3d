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
   cout << "volu : index " << index << " | fill " << fill << ", " << " | isBorder " << isBorder <<
   target << ", "<< error << " | surfind " << surfind.size();
   for (i=0; unsigned_int(i)<surfind.size();i++){
      cout << "-" << surfind[i];
   }
   cout << endl;
}

void surf::disp() const{
   int i;
   cout << "surf : index " << index << " | fill " << fill << " | isBorder " << isBorder << ", " << 
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
   bool needRePrep;

   verts.PrepareForUse();
   edges.PrepareForUse();
   surfs.PrepareForUse();
   volus.PrepareForUse();
   // Additional mesh preparation steps

   needRePrep=this->OrderEdges();
   if (needRePrep){
      surfs.PrepareForUse();
   }

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
      edgeIndOrig=edgeind;
      edgeSub=meshin->edges.find_list(edgeind);
      edge2Vert=ConcatenateVectorField(meshin->edges,&edge::vertind,edgeSub);

      HashVector(edge2Vert, vert2Edge);

      vertCurr=edge2Vert[0];
      edgeCurr=edgeind[0];
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
      #endif // SAFe_ACCESS
         jj=edgeIndOrig[(range.first->second)/2]==edgeCurr;

         it=range.first;
         if (jj){++it;}

         edgeCurr=edgeIndOrig[(it->second)/2];
         jj=edge2Vert[((it->second)/2)*2]==vertCurr; // Warnign ((x/2)*2) not necessarily equal
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
   // Update border status of volus
   for(ii=0;ii<volus.size();++ii){
      jj=0;
      nT=volus(ii)->surfind.size();
      while(jj<nT && !volus(ii)->isBorder){
         volus[ii].isBorder=surfs(surfs.find(volus[ii].surfind[jj]))->isBorder;
         jj++;
      }
   }
   // Update border status of edges
   for(ii=0;ii<edges.size();++ii){
      jj=0;
      nT=edges(ii)->surfind.size();
      while(jj<nT && !edges(ii)->isBorder){
         edges[ii].isBorder=surfs(surfs.find(edges[ii].surfind[jj]))->isBorder;
         jj++;
      }
   }
   // Update border status of edges
   for(ii=0;ii<verts.size();++ii){
      jj=0;
      nT=verts(ii)->edgeind.size();
      while(jj<nT && !verts(ii)->isBorder){
         verts[ii].isBorder=surfs(edges.find(verts[ii].edgeind[jj]))->isBorder;
         jj++;
      }
   }


}

