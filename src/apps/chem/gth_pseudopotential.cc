
#include <chem/gth_pseudopotential.h>


namespace madness {


    double get_charge_from_file(const std::string filename, unsigned int atype) {
        bool debug = true;
       
        TiXmlDocument doc(filename);
        if (!doc.LoadFile()) {
            MADNESS_EXCEPTION("Failed to load GTH pseudopotential file", 0);
        }

        bool success = false;
        for (TiXmlElement* node=doc.FirstChildElement(); node && !success; node=node->NextSiblingElement()) {
            if (strcmp(node->Value(),"name") == 0) {
                std::string name = node->GetText();
                if (debug) std::cout << "Loading pseudopotential file " << name << std::endl;
            }
            else if (strcmp(node->Value(), "atom") == 0) {
                const char* symbol = node->Attribute("symbol");
                unsigned int atn = symbol_to_atomic_number(symbol);
                if (atype == atn) {
                    success = true;
                    if (debug) std::cout << "  found atomic pseudopotential " << symbol << std::endl;
                    int lmax = -1;
                    node->Attribute("lmax", &lmax);
                    // local part
                    TiXmlElement* xmlVLocal = node->FirstChildElement();
                    double zeff = 0.0; xmlVLocal->Attribute("Zeff", &zeff);
                    return zeff;
                }
            }
        }
        if (!success){
            MADNESS_EXCEPTION("Failed to find element in GTH pseudopotential file", 0);
        }
        
        return 0.;
    }

    const double ProjRLMFunctor::gamma_data[17] = {1.0, 0.0, 1.0/2.0, 0.0, 3.0/4.0, 0.0, 15.0/8.0, 0.0, 105.0/16.0, 0.0, 
        945.0/32.0, 0.0, 10395.0/64.0, 0.0, 135135.0/128.0, 0.0, 2027025.0/256.0};


}
