// diag_times_right.cpp
// MATLAB mex function that multiplies a stack of matricies H by a diagnoal matrix vec from the left. 
// The diagonal matrix vec is names so because only the diagonals are saved of the diagonal matrix.
// 
// USAGE:
// M_added = add_diag(M, v);
#include "mex.hpp"
#include "mexAdapter.hpp"
#include <iostream>

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
        TypedArray<double> H = std::move(inputs[0]);
        TypedArray<double> vec = std::move(inputs[1]);
        for(auto i = 0; i < H.getDimensions()[0]; i++){
            for(auto j = 0; j < H.getDimensions()[1]; j++){
                for(auto k = 0; k < H.getDimensions()[2]; k++){
                    H[i][j][k] *= vec[j][k];
                }
            }
        }
        outputs[0] = H;
    }

    //Argument checks for the inputs. I don't really understand why MATLAB says we should pass in the 
    //output arguments.
    void checkArguments(ArgumentList outputs, ArgumentList inputs) {
        if (inputs[0].getType() != ArrayType::DOUBLE)
        {
            matlabPtr->feval(u"error", 0, 
                    std::vector<Array>({ factory.createScalar("M must be of type double") }));
        }
        if (inputs[0].getDimensions()[0] !=  inputs[0].getDimensions()[1])
        {
            matlabPtr->feval(u"error", 0, 
                    std::vector<Array>({ factory.createScalar("Must be square matrix") }));
        }
        if (inputs[1].getType() != ArrayType::DOUBLE)
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
