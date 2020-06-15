#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h> 
#include <sys/types.h>
#include <algorithm>
#include <vector>

using namespace std;

class Answer
{
public:
	static bool exists(int ints[], int size, int k)
	{
		vector<int> vec;
		for (int i = 0; i < size; ++i) for (int* ptrInt = ints; ptrInt < ints + size; ptrInt++) {
			vec.push_back(move(ints[i]));
		}
		if (binary_search(vec.begin(), vec.end(), k)) {
			return true;
		}
		else {
			return false;
		}

	}
};
int main() {
	int ints[] = { 4, 6 ,7 ,8 };
	cout << Answer::exists(ints, 4, 4) << endl;
	cout << Answer::exists(ints, 4, 10) << endl;


}