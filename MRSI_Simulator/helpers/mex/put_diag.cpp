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
        for(int j = 0; j < H.getDimensions()[2]; j ++){
            for(int i = 0; i < vec.getDimensions()[0]; i++){
                H[i][i][j] = vec[i][j];
            }
        }
        outputs[0] = H;
    }

    //Argument checks for the inputs. I don't really understand why MATLAB says we should pass in the 
    //output arguments.
    void checkArguments(ArgumentList outputs, ArgumentList inputs) {
        if (inputs[0].getType() != ArrayType::COMPLEX_SINGLE)
        {
            matlabPtr->feval(u"error", 0, 
                    std::vector<Array>({ factory.createScalar("M must be of type double") }));
        }
        if (inputs[0].getDimensions()[0] !=  inputs[0].getDimensions()[1])
        {
            matlabPtr->feval(u"error", 0, 
                    std::vector<Array>({ factory.createScalar("Must be square matrix") }));
        }
        if (inputs[1].getType() != ArrayType::COMPLEX_SINGLE)
        {
            matlabPtr->feval(u"error", 0, 
                    std::vector<Array>({ factory.createScalar("Vec must be of type double") }));
        }
        if (inputs[1].getDimensions()[0] != inputs[0].getDimensions()[0])
        {
            matlabPtr->feval(u"error", 0, 
                    std::vector<Array>({ factory.createScalar("Vec and H must have the same first dimension size") }));
        }
        if (inputs[1].getDimensions()[1] != inputs[0].getDimensions()[2])
        {
            matlabPtr->feval(u"error", 0, 
                    std::vector<Array>({ factory.createScalar("Vec second dimension and H third dimensions need to be the same sizes") }));
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
