#include "rsvsclasses2json.hpp"

void rsvsjson::flatupdate(json& jfin, json& jnew,
	bool isFlatFin, bool isFlatNew){
	/*
	Allows recursing updae into sub fields
	*/
	// std::cout << "file read " << std::endl;
	if(!isFlatNew){
		jnew = jnew.flatten();
	}
	if(!isFlatFin){
		jfin = jfin.flatten();
	}
		// std::cout << "j unflattened " << std::endl;
		// std::cout << jnew.dump(1) << std::endl;
		// std::cout << jfin.dump(1) << std::endl;
	// Insert values read into the parameter structure
	jfin.update(jnew);
		// std::cout << "jfin updated " << std::endl;
		// std::cout << jfin.dump(1) << std::endl;
	jfin = jfin.unflatten();
	jnew = jnew.unflatten();
}