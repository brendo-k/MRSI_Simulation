// exp_gpu.cpp
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
    }

    //Argument checks for the inputs. I don't really understand why MATLAB says we should pass in the 
    //output arguments.
    void checkArguments(ArgumentList outputs, ArgumentList inputs) {
        if (inputs[0].getType() != ArrayType::DOUBLE)
        {
            matlabPtr->feval(u"error", 0, 
                    std::vector<Array>({ factory.createScalar("M must be of type double") }));
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
