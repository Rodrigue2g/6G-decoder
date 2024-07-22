% This files is used to compile all mex functions

% Start the environment
mex -setup;

% % Floating-point
mex -output decoding\mexw64\sc_ldpc_layered_nms_float_decoding_mex decoding\sc_ldpc_layered_nms_float_decoding_mex.c;
mex -output decoding\mexw64\sc_ldpc_flooding_nms_float_decoding_mex decoding\sc_ldpc_flooding_nms_float_decoding_mex.c;