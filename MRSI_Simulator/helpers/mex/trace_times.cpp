
// put_diag.cpp
// MATLAB mex function that puts values along diagonals. 
// 
// USAGE:
// M_put = put_diag(M, v);
#include "mex.hpp"
#include "mexAdapter.hpp"
#include <iostream>
#include <utility>
#include <vector>

using namespace matlab::data;
using matlab::mex::ArgumentList;

class MexFunction : public matlab::mex::Function {
    std::shared_ptr<matlab::engine::MATLABEngine> matlabPtr = getEngine();
    // Factory to create MATLAB data arrays
    ArrayFactory factory;
    // Create an output stream
    std::ostringstream stream;
    

    public: 
    //Entry point for mex function. This is what is called when calling from matlab
    void operator()(ArgumentList outputs, ArgumentList inputs){
        checkArguments(outputs, inputs);
        TypedArray<std::complex<float>> H = std::move(inputs[0]);
        TypedArray<std::complex<float>> vec = std::move(inputs[1]);
        int sum = 0;
        for (int i = 0; int i < H.getDimensions()[2]; i++) {
            for (int j = 0; j < H.getDimensions()[1]; j++){
                sum += 
            }
        }
        
    }

    //Argument checks for the inputs. I don't really understand why MATLAB says we should pass in the 
    //output arguments.
    void checkArguments(ArgumentList outputs, ArgumentList inputs) {
        if (inputs[0].getType() != ArrayType::COMPLEX_SINGLE)
        {
            matlabPtr->feval(u"error", 0, 
                    std::vector<Array>({ factory.createScalar("Spins must be of type single") }));
        }
        if (inputs[0].getDimensions()[0] !=  inputs[0].getDimensions()[1])
        {
            matlabPtr->feval(u"error", 0, 
                    std::vector<Array>({ factory.createScalar("Must be square matrix") }));
        }
        if (inputs[1].getType() != ArrayType::COMPLEX_SINGLE)
        {
            matlabPtr->feval(u"error", 0, 
                    std::vector<Array>({ factory.createScalar("Matrix must be of type single") }));
        }
        if (inputs[0].getDimensions()[0] !=  inputs[0].getDimensions()[1])
        {
            matlabPtr->feval(u"error", 0, 
                    std::vector<Array>({ factory.createScalar("Marix must be square") }));
        }
        if (inputs[0].getDimensions()[1] != inputs[1].getDimensions()[0]) {
            matlabPtr->feval(u"error", 0, 
                    std::vector<Array>({ factory.createScalar("Matrix dimensions must agree")}));
        }
        

    }
    void displayOnMATLAB(std::ostringstream& stream) {
        // Pass stream content to MATLAB fprintf function
        matlabPtr->feval(u"fprintf", 0,
                std::vector<Array>({ factory.createScalar(stream.str()) }));
        // Clear stream buffer
        stream.str("");
    }
};
