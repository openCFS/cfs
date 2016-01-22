
#include "coeff.h"
#include <iostream>
#include <fstream>
using namespace std;


void write_matrix() {
	ofstream output("matrix.dat");
	for(int i=0;i<64;i++) {
		for(int j=0;j<64;j++){
			output<<A[i][j]<<" ";
		}
		output<<endl;
	}
}

int main(void) {
write_matrix();
return 0;
}
