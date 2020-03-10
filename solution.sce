// Task 1 {{{1 =================================================================

// Custom convolution function {{{2 ============================================

// This function uses the multiplicative property of FFT to compute convolution
// with the speed of an FFT (which is NlogN instead of N^2).
function result = custom_convolve(x, y)
  zero_padded_x = cat(2, x, zeros(1, length(y) - 1))
  zero_padded_y = cat(2, y, zeros(1, length(x) - 1))
  result = ifft(fft(zero_padded_x) .* fft(zero_padded_y))
endfunction // custom_convolve

// Loading the IRC {{{2 ========================================================

irc = loadwave('irc.wav')

// IRC was recorded at 44100 Hz, and the default is 22050 Hz, so taking every
// 2-nd value
irc = irc(1 : 2 : length(irc))

// Convolving the IRC with the sound using the custom function {{{2 ============

// File naming convention for the following: `proc_` for the reference outputs,
// `output_` for the custom outputs.

violin = loadwave('violin.wav')
custom_new_violin = custom_convolve(irc, violin)
wavwrite(custom_new_violin, 'output_violin.wav')

drums = loadwave('drums.wav')
custom_new_drums = custom_convolve(irc, drums)
wavwrite(custom_new_drums, 'output_drums.wav')

voice = loadwave('voice.wav')
custom_new_voice = custom_convolve(irc, voice)
wavwrite(custom_new_voice, 'output_voice.wav')

// This one has two channels.
// Let's take only one of them for the sake of simplicity.
speech = loadwave('speech.wav')
speech = speech(1,:)
custom_new_speech = custom_convolve(irc, speech)
wavwrite(custom_new_speech, 'output_speech.wav')

// Convolving the IRC with the sound using the built-in function {{{2 ==========

reference_new_violin = convol(irc, violin)
wavwrite(reference_new_violin, 'proc_violin.wav')

reference_new_drums = convol(irc, drums)
wavwrite(reference_new_drums, 'proc_drums.wav')

reference_new_voice = convol(irc, voice)
wavwrite(reference_new_voice, 'proc_voice.wav')

reference_new_speech = convol(irc, speech)
wavwrite(reference_new_speech, 'proc_speech.wav')

// Task 2. IIR {{{1 ============================================================

// Function implementation {{{2 ================================================

// Constructor
function result = init_iir(recursive_coeffs, source_coeffs)
  result = struct()

  // Recent outputs of the IIR (y in the formula)
  result.recent_results = zeros(recursive_coeffs)

  // Recent inputs to the IIR (x in the formula)
  result.recent_inputs = zeros(source_coeffs)

  // Coefficients for the recursive part (a in the formula)
  result.recursive_coeffs = flipdim(recursive_coeffs, 2)

  // Coefficients for the source part (b in the formula)
  result.source_coeffs = flipdim(source_coeffs, 2)
endfunction // init_iir

// Advance IIR by one value
// (modifying the internal state, and returning a new IIR)
function result = advance_iir(iir, input_value)
  result = iir

  // Add an input into the recent_inputs
  shifted_recent_inputs = result.recent_inputs(2:length(result.recent_inputs))
  result.recent_inputs = cat(2, shifted_recent_inputs, [input_value])

  // Calculate the output value ...
  output_value  ..
    = sum(result.recursive_coeffs .* result.recent_results)  ..
    + sum(result.source_coeffs .* result.recent_inputs)

  // ... and append it to the struct
  shifted_recent_results = result.recent_results(2:length(result.recent_results))
  result.recent_results = cat(2, shifted_recent_results, [output_value])
endfunction // advance_iir

// Compute the IIR over the input signal (destroying the IIR)
function result = iir_of(iir, input_values)
  output_values = []
  for i = [1 : length(input_values)]
    iir = advance_iir(iir, input_values(i))
    output = iir.recent_results(length(iir.recent_results))
    output_values = cat(2, output_values, [output])

    if modulo(i, 1000) == 0 then
      mprintf('%i out of %i\n', i, length(input_values))
    end
  end

  result = output_values
endfunction

// Using the IIR to lowpass the sound {{{2 =====================================

violin_sound = loadwave('Violin_Viola_Cello_Bass.wav')

highpass = init_iir(  ..
  [-0.3769782747249014, -0.19680764477614976],  ..
  [0.40495734254626874, -0.8099146850925375, 0.4049573425462687]  ..
)
high_violin_sound = iir_of(highpass, violin_sound)
wavwrite(high_violin_sound, 'output_high.wav')

lowpass = init_iir(  ..
  [1.9733442497812987, -0.9736948719763],  ..
  [0.00008765554875401547, 0.00017531109750803094, 0.00008765554875401547]  ..
)
low_violin_sound = iir_of(lowpass, violin_sound)
wavwrite(low_violin_sound, 'output_low.wav')

// vim: set fdm=marker :
