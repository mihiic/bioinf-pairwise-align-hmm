#include <iostream>
#include <fstream>
#include "alignment_model.h"

int main(int argc, char** argv) {
    AlignmentModel model(60);
    if (argc != 3) {
        cout << "You have to specify exactly 2 input file." << endl;
        return 1;
    }

    model.run(argv[1], argv[2]);
    return 0;
}