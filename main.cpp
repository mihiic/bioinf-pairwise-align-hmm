#include <iostream>
#include <fstream>
#include "alignment_model.h"

int main(int argc, char** argv) {
    AlignmentModel model(60);
    model.run(argv[1], argv[2]);
    return 0;
}