#include <vector>
#include <iostream>

using namespace std;

int main(int argc, const char* argv[]) {
	vector < vector <int> > psi;
	vector < vector <float> > delta;
	vector < string* > nodeList;
	
	size_t tmp = 144;
	
	for (int i=0; i<500; i++) {
		nodeList.push_back(new string("This is a test"));
		psi.push_back(vector <int> (tmp,0));
		cout << " Psi cap: " << psi.capacity() << endl;
		delta.push_back(vector <float> (tmp,0));
		cout << " Delta cap: " << delta.capacity() << endl;
		for (int j=0; j<10; j++) {
			psi[i][j]=j;
			delta[i][j]=j*5.5;
			//psi[i].push_back(j);
			//delta[i].push_back(j*5.5);
		}
	}
	for (int i=0; i<500; i++) {
		for (int j=0; j<10; j++) {
			cout << psi[i][j] << " ";
		}
		cout << endl;
		for (int j=0; j<10; j++) {
			cout << delta[i][j] << " ";
		}
		cout << endl;
	}
}
