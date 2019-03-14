#include "rsvsjson.hpp"

void rsvsjson::flatupdate(rsvsjson::json& jfin, rsvsjson::json& jnew,
	bool isFlatFin, bool isFlatNew){
	/*
	Allows recursing update into sub fields
	*/
	if(!isFlatNew){
		jnew = jnew.flatten();
	}
	if(!isFlatFin){
		jfin = jfin.flatten();
	}

	// Insert values read into the parameter structure
	jfin.update(jnew);

	jfin = jfin.unflatten();
	jnew = jnew.unflatten();
}