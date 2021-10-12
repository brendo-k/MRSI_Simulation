#include "mex.hpp"
#include "mexAdapter.hpp"
#include <iostream>
#include <utility>

using namespace matlab::data;
using matlab::mex::ArgumentList;

class MexFunction : public matlab::mex::Function {
    std::shared_ptr<matlab::engine::MATLABEngine> matlabPtr = getEngine();
    // Factory to create MATLAB data arrays
    ArrayFactory factory;
    // Create an output stream
    std::ostringstream stream;


public: 
    void operator()(ArgumentList outputs, ArgumentList inputs){
       checkArguments(outputs, inputs); 
       TypedArray<Struct> mets = std::move(inputs[0]);
       TypedArray<std::complex<float>> gradients = std::move(inputs[1]);
       TypedArray<double> spins = matlabPtr -> getProperty(mets, u"spins");

       float gamma = 42577000;
       calculate_ppm(gradients, spins);

    }
    void checkArguments(ArgumentList outputs, ArgumentList inputs) {
        if (inputs[0].getType() != ArrayType::STRUCT)
        {
            matlabPtr->feval(u"error", 0, 
                std::vector<Array>({ factory.createScalar("Phantom input must be a struct") }));
        }

        if (inputs[1].getType() != ArrayType::DOUBLE || inputs[1].getDimensions().size() != 2)
        {
            matlabPtr->feval(u"error", 0, 
                std::vector<Array>({ factory.createScalar("Gradient input must be double array") }));
        }
        

        if (inputs[2].getType() != ArrayType::DOUBLE || inputs[2].getNumberOfElements() != 1)
        {
            matlabPtr->feval(u"error", 0, 
                std::vector<Array>({ factory.createScalar("time must be a scalar") }));
        }

        if (outputs.size() > 6) {
            matlabPtr->feval(u"error", 0, 
                std::vector<Array>({ factory.createScalar("Not enough inputs") }));
        }
    }
    void displayOnMATLAB(std::ostringstream& stream) {
        // Pass stream content to MATLAB fprintf function
        matlabPtr->feval(u"fprintf", 0,
            std::vector<Array>({ factory.createScalar(stream.str()) }));
        // Clear stream buffer
        stream.str("");
    }

    void calculate_ppm(TypedArray<std::complex<float>>& gradients, 
                                                    TypedArray<double> spins)
    {
        TypedArray<std::complex<float>> temp = TypedArray<std::complex<float>>(gradients);
        Array result = matlabPtr -> feval(u"bsxfun", "@times", gradients, factory.createScalar(spins[0]/1000000));
        Array result2 = matlabPtr -> feval(u"bsxfun", "@plus", temp, factory.createScalar(-1));
        Array result3 = matlabPtr -> feval(u"plus", result, result2);
        Array result4 = matlabPtr -> feval(u"bsxfun", "@times", result3, factory.createScalar(1000000));
        gradients = std::move(result4);
    }
};
