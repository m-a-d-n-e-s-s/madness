//
// Created by Florian Bischoff/chatgpt on 9/7/23.
//

#include<madness/misc/misc.h>
#include<madness/misc/ran.h>

namespace madness {

std::string unique_fileid() {
    std::string uniqueFileName;

    // Check if the PBS_JOBID environment variable is set
    const char* pbsId = std::getenv("PBS_JOBID");
    if (pbsId != nullptr) {
        uniqueFileName =std::string(pbsId);
    } else {
        // If PBS_ID is not available, generate a random number
        int randomNumber = RandomValue<int>();
        uniqueFileName = std::to_string(randomNumber);
    }

    return uniqueFileName;
}
}
