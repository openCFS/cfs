#include <string>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <iomanip>

// #include <muParser.h>

#include "params.hh"
#include "filereader.hh"


namespace CoupledField
{

    FileReader::FileReader(const std::string& name,
                           const UInt dim,
                           const UInt numFiles) :
        numPartitions_(0), 
        dim_(dim),
        numResults_(1),
        numFiles_(numFiles)
    {
        name_ = name;
        basename_ = "./";
        basename_+= name_;
        basename_+= "_data/";
        basename_+= name_;
        basename_+= "_";
    }

    FileReader::~FileReader()
    {
    }

    double FileReader::GetTimeStep()
    {
        double ts;
/*        
        mu::Parser parser;

        parser.SetExpr(params->GetTimeStep());
        try
        {
            ts = parser.Eval();
        }
        catch (mu::Parser::exception_type &e)
        {
            std::cerr << e.GetMsg() << std::endl;
            exit(1);
        }

        if(ts < 0)
        {
            std::cerr << "Negative timestep!"
                      << "You must specify timestep with -p or --timestep parameters."
                      << std::endl;
            exit(1);
        }
*/
        return ts;
    }
      
} // end of namespace
