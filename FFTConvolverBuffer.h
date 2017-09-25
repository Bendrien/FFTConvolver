//
// Created by anyere on 18.09.17.
//

#ifndef FFTCONVOLVER_FFTBUFFER_H
#define FFTCONVOLVER_FFTBUFFER_H

#include "AudioFFT.h"
#include "Utilities.h"

#include <vector>

namespace fftconvolver {

class FFTConvolverBuffer {
public:
    FFTConvolverBuffer();
    virtual ~FFTConvolverBuffer();

    /**
    * @brief Initializes the convolver
    * @param blockSize Block size internally used by the convolver (partition size)
    * @param ir The impulse response
    * @param irLen Length of the impulse response
    * @return true: Success - false: Failed
    */
    bool init(size_t blockSize, const Sample* ir, size_t irLen);

    /**
    * @brief Convolves the the given input samples and immediately outputs the result
    * @param input The input samples
    * @param output The convolution result
    * @param len Number of input/output samples
    */
    void process(const Sample* input, Sample* output, size_t len);

    /**
    * @brief Resets the convolver and discards the set impulse response
    */
    void reset();

private:
    /// Block size internally used by the convolver (partition size, is a power of 2)
    size_t _blockSize;

    /// Size of one segment in samples (block size times 2)
    size_t _segSize;

    /// Count of all segments (fitting inside the impulse response)
    size_t _segCount;

    /// Size of real/imagenary part of the segments in the frequency domain
    size_t _fftComplexSize;

    /// Vector of all buffered input segments in the frequency domain
    std::vector<SplitComplex*> _segments;

    /// Vector of all impulse response segments in the frequency domain
    std::vector<SplitComplex*> _segmentsIR;

    /// Sample buffer, holding the converted input into the frequency domain (sized by segment size)
    SampleBuffer _fftBuffer;

    /// AudioFFT handle
    audiofft::AudioFFT _fft;

    /// Buffer for convolution proposes (sized by segment size)
    SplitComplex _preMultiplied;

    /// Buffer for convolution proposes (sized by segment size)
    SplitComplex _conv;

    /// Sample buffer, holding the overlap (sized by buffer size)
    SampleBuffer _overlap;

    /// Index of the current segment
    size_t _current;

    /// Sample buffer, holding the converted input into the frequency domain (sized by buffer size)
    SampleBuffer _inputBuffer;

    /// Position of the processed input in samples
    size_t _inputBufferFill;

    // Prevent uncontrolled usage
    FFTConvolverBuffer(const FFTConvolverBuffer&);
    FFTConvolverBuffer& operator=(const FFTConvolverBuffer&);
};

} // End of namespace fftconvolver

#endif //FFTCONVOLVER_FFTBUFFER_H
