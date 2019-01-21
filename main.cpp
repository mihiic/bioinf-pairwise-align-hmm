#include <iostream>
#include <fstream>
#include <string>
#include "alignment_model.h"

int main(int argc, char** argv) {
    if (argc == 7) {
        AlignmentModel model((unsigned long)stoi(argv[6]));
        model.run(argv[1], argv[2], argv[3], argv[4], argv[5]); //
    } else {
        AlignmentModel model(4096);
        model.run(argv[1], argv[2], argv[3], argv[4], argv[5]); //
    }
    return 0;
}