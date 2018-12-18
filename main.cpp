#include <iostream>
#include <fstream>
#include "alignment_model.h"

int main(int argc, char** argv) {
    AlignmentModel model(4096);
    model.run(argv[1], argv[2], argv[3], argv[4], argv[5]); //
    return 0;
}