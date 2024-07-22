library ieee;
use ieee.std_logic_1164.all;
use ieee.numeric_std.all;
use ieee.math_real.all;



package constants is

    -- Size constants
        -- MxN base graph matrix (BG1 BG2)
        constant N_BG : integer := 80;          
        constant M_BG : integer := 8;

        -- Z is 16
        constant Z : integer := 16;

        -- (M*Z)x(N*Z) full LDPC matrix (H1 H2)
        constant N : integer := Z*N_BG;                     -- this is also the size of the codewords used in the decoding process
        constant M : integer := Z*M_BG;

        -- total number of codewords
        constant L : integer := 50; -- size of SC-LDPC code
        
        -- (M*Z*L)x(N*Z*L) full SC-LDPC matrix 
        constant N_FULL : integer := N*L;
        constant M_FULL : integer := M*L;

        -- number of bits needed to represent the constants
        constant N_BG_bits : integer := integer(ceil(log2(real(N_BG))));
        constant M_BG_bits : integer := integer(ceil(log2(real(M_BG))));
        constant Z_bits : integer := integer(ceil(log2(real(Z))))+1;  -- + 1 because we also need to represent -1
        constant N_bits : integer := integer(ceil(log2(real(N))));
        constant M_bits : integer := integer(ceil(log2(real(M))));
        constant L_bits : integer := integer(ceil(log2(real(L))));
        constant N_FULL_bits : integer := integer(ceil(log2(real(N_FULL))));
        constant M_FULL_bits : integer := integer(ceil(log2(real(M_FULL))));

        -- Total number of bits in the base graph (that we need to fit in the FPGA)
        constant BG_SIZE : integer := N_BG*M_BG*Z_bits;
        


    -- Constants for fixed point representation

        -- the fixed point representation is Q3.5 (8 bits)
        constant FXP_WIDTH : integer := 8; --with sign bit
        constant num_integer_bits : natural := 3; -- Number of integer bits (including sign bit)
        constant num_fractional_bits : natural := 5; -- Number of fractional bits
        constant scale_factor : real := 2.0**num_fractional_bits; -- Scaling factor

        constant METRIC_WIDTH : integer := 17; --with sign bit 
        constant num_metric_integer_bits : natural := 12; -- Number of integer bits (including sign bit)
        constant num_metric_fractional_bits : natural := 5; -- Number of fractional bits
        constant metric_scale_factor : real := 2.0**num_metric_fractional_bits; -- Scaling factor

    -- Data types
        type t_BG_matrix is array(0 to M_BG-1, 0 to N_BG-1) of signed(Z_bits-1 downto 0);
        type t_y_array is array(0 to N-1) of signed(FXP_WIDTH-1 downto 0);
        type t_y_array_half is array(0 to N/2-1) of signed(FXP_WIDTH-1 downto 0);
        type t_xhat_array is array(0 to N-1) of std_logic;
        type t_xhat_array_half is array(0 to N/2-1) of std_logic;
        subtype t_metric is signed(METRIC_WIDTH-1 downto 0);


    -- Constants for the LDPC decoder
        constant MAX_ITERATIONS     : integer := 100;    -- Maximum number of iterations of the LDPC decoder
        constant BETA               : t_metric := "11111111111110000";  -- -0.5 in fixed point representation Q3.5

        constant MAX_ITERATIONS_bits    : integer := integer(ceil(log2(real(MAX_ITERATIONS))));


end package constants;
