library ieee;
use ieee.std_logic_1164.all;
use ieee.numeric_std.all;

library work;
use work.constants.all;

-- ! testbench
--use ieee.std_logic_textio.all;
entity MetricCalc is
    port
    (
        CLKxCI : in std_logic;
        RSTxRI : in std_logic;

        xhatxDI : in t_xhat_array;
        yxDI    : in t_y_array;

        startxSI : in std_logic;

        validxSO  : out std_logic;
        parityxSO : out std_logic;
        valuexDO  : out t_metric;

        flipvalidxSO : out std_logic;
        flipindexxDO : out unsigned(N_bits - 1 downto 0);

        SingleFlipindexxDO : out unsigned(N_bits - 1 downto 0)
    );
end entity MetricCalc;
architecture rtl of MetricCalc is

    -- the basegraph
    signal BG : t_BG_matrix;

    --IDLE waits for start signal, STEP_ONE takes 128 cycles, STEP_TWO takes 1280 cycles
    type t_state is (IDLE, STEP_ONE, STEP_TWO, END_DEC);

    -- holds the result of the first step
    signal VecOnexDN, VecOnexDP : std_logic_vector(M - 1 downto 0);

    -- used for the objective function calculation
    signal sum_paritiesxDN, sum_paritiesxDP       : signed(M_bits + 1 downto 0); -- M_bits+2 because of the sign bit and the fact we can reach 128
    signal sum_correlationxDN, sum_correlationxDP : t_metric;

    -- counter that goes up to the length of the codeword (1280)
    signal counterxDN, counterxDP : unsigned(N_bits - 1 downto 0);

    -- control signals
    signal validxDN, validxDP                       : std_logic;
    signal parityxDN, parityxDP                     : std_logic;
    signal argmin_inverfuncxDN, argmin_inverfuncxDP : t_metric;
    signal argmin_indexxDN, argmin_indexxDP         : unsigned(N_bits - 1 downto 0);
    signal statexDN, statexDP                       : t_state;
    signal previndexxDN, previndexxDP               : unsigned(N_bits - 1 downto 0);
    signal term1xDN, term1xDP                       : unsigned(M_BG_bits-1 downto 0);
    signal term2xDN, term2xDP                       : unsigned(M_BG_bits-1 downto 0);
    signal term3xDN, term3xDP                       : t_metric;
    signal flipvalidxS: std_logic;

begin

    BG <= (
        (to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(13, Z_bits), to_signed(-1, Z_bits), to_signed(14, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(2, Z_bits), to_signed(5, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(9, Z_bits), to_signed(3, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(4, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(3, Z_bits), to_signed(-1, Z_bits), to_signed(13, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(13, Z_bits), to_signed(-1, Z_bits), to_signed(11, Z_bits), to_signed(7, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(9, Z_bits), to_signed(14, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(9, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(6, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(14, Z_bits), to_signed(2, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(0, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits)),
        (to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(12, Z_bits), to_signed(-1, Z_bits), to_signed(6, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(1, Z_bits), to_signed(8, Z_bits), to_signed(-1, Z_bits), to_signed(14, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(13, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(10, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(4, Z_bits), to_signed(-1, Z_bits), to_signed(8, Z_bits), to_signed(4, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(5, Z_bits), to_signed(-1, Z_bits), to_signed(2, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(11, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(6, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(10, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(12, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(4, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(12, Z_bits), to_signed(-1, Z_bits), to_signed(0, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits)),
        (to_signed(14, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(9, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(0, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(2, Z_bits), to_signed(-1, Z_bits), to_signed(10, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(0, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(15, Z_bits), to_signed(13, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(8, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(3, Z_bits), to_signed(5, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(14, Z_bits), to_signed(9, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(2, Z_bits), to_signed(-1, Z_bits), to_signed(0, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(8, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(5, Z_bits), to_signed(0, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(3, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(0, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits)),
        (to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(3, Z_bits), to_signed(-1, Z_bits), to_signed(10, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(11, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(12, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(9, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(10, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(0, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(3, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(11, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(1, Z_bits), to_signed(1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(5, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(11, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(2, Z_bits), to_signed(15, Z_bits), to_signed(-1, Z_bits), to_signed(8, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(12, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(2, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(3, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(0, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits)),
        (to_signed(5, Z_bits),  to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(3, Z_bits),  to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(11, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(7, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(5, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(5, Z_bits), to_signed(-1, Z_bits), to_signed(5, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(9, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(2, Z_bits), to_signed(1, Z_bits), to_signed(1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(13, Z_bits), to_signed(13, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(15, Z_bits), to_signed(9, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(6, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(3, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(1, Z_bits), to_signed(13, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(0, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits)),
        (to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(5, Z_bits),  to_signed(1, Z_bits),  to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(6, Z_bits), to_signed(-1, Z_bits), to_signed(3, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(14, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(13, Z_bits), to_signed(-1, Z_bits), to_signed(15, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(14, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(2, Z_bits), to_signed(13, Z_bits), to_signed(12, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(5, Z_bits), to_signed(1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(7, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(11, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(0, Z_bits), to_signed(6, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(13, Z_bits), to_signed(-1, Z_bits), to_signed(7, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(0, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits)),
        (to_signed(-1, Z_bits), to_signed(11, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(6, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(14, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(12, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(13, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(10, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(10, Z_bits), to_signed(-1, Z_bits), to_signed(4, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(11, Z_bits), to_signed(9, Z_bits), to_signed(-1, Z_bits), to_signed(13, Z_bits), to_signed(-1, Z_bits), to_signed(11, Z_bits), to_signed(-1, Z_bits), to_signed(1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(7, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(1, Z_bits), to_signed(1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(4, Z_bits), to_signed(-1, Z_bits), to_signed(5, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(0, Z_bits), to_signed(-1, Z_bits)),
        (to_signed(-1, Z_bits), to_signed(7, Z_bits),  to_signed(10, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(7, Z_bits), to_signed(2, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(0, Z_bits), to_signed(11, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(15, Z_bits), to_signed(-1, Z_bits), to_signed(11, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(9, Z_bits), to_signed(8, Z_bits), to_signed(11, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(0, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(13, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(10, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(14, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(9, Z_bits), to_signed(3, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(12, Z_bits), to_signed(-1, Z_bits), to_signed(11, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(0, Z_bits))
        );

    p_dff : process (CLKxCI, RSTxRI)
    begin
        if RSTxRI = '1' then
            validxDP            <= '0';
            parityxDP           <= '0';
            counterxDP          <= (others => '0');
            sum_paritiesxDP     <= (others => '0');
            sum_correlationxDP  <= (others => '0');
            VecOnexDP           <= (others => '0');
            argmin_inverfuncxDP <= (others => '0');
            argmin_indexxDP     <= (others => '0');
            statexDP            <= IDLE;
            previndexxDP        <= (others => '0');
            term1xDP            <= (others => '0');
            term2xDP            <= (others => '0');
            term3xDP            <= (others => '0');



        elsif rising_edge(CLKxCI) then
            validxDP            <= validxDN;
            parityxDP           <= parityxDN;
            counterxDP          <= counterxDN;
            sum_paritiesxDP     <= sum_paritiesxDN;
            sum_correlationxDP  <= sum_correlationxDN;
            VecOnexDP           <= VecOnexDN;
            argmin_inverfuncxDP <= argmin_inverfuncxDN;
            argmin_indexxDP     <= argmin_indexxDN;
            statexDP            <= statexDN;
            previndexxDP        <= previndexxDN;
            term1xDP            <= term1xDN;
            term2xDP            <= term2xDN;
            term3xDP            <= term3xDN;


        end if;
    end process p_dff;
    p_main : process (all)
        variable line_vec1 : std_logic_vector(0 to N_BG - 1);
        variable vec1_xor  : std_logic;
        variable xtimesy   : t_metric;
        variable line_vec2 : std_logic_vector(0 to M_BG - 1);
        variable line_vec3 : std_logic_vector(0 to M_BG - 1);
        variable sum_vec2  : unsigned(M_BG_bits-1 downto 0);
        variable sum_vec3  : unsigned(M_BG_bits-1 downto 0);

        variable a1 : signed(Z_bits-1 downto 0);
        variable a2 : unsigned(Z_bits-2 downto 0);
        variable b1 : unsigned(Z_bits-2 downto 0);
        variable y_resized_shifted : t_metric;

    begin

        -- default values
        validxDN            <= '0';
        parityxDN           <= parityxDP;
        counterxDN          <= counterxDP;
        sum_paritiesxDN     <= sum_paritiesxDP;
        sum_correlationxDN  <= sum_correlationxDP;
        VecOnexDN           <= VecOnexDP;
        statexDN            <= statexDP;
        argmin_inverfuncxDN <= argmin_inverfuncxDP;
        argmin_indexxDN     <= argmin_indexxDP;
        previndexxDN        <= previndexxDP;
        term1xDN            <= term1xDP;
        term2xDN            <= term2xDP;
        term3xDN            <= term3xDP;

        -- the second pipeline stage of the inversion function where the values are summed and compared (to reduce critical path)
        xtimesy := term3xDP + signed(shift_left(resize(term2xDP, xtimesy'length), num_metric_fractional_bits)) - signed(shift_left(resize(term1xDP, xtimesy'length), num_metric_fractional_bits));
        if  xtimesy < argmin_inverfuncxDP then
            argmin_inverfuncxDN <=  xtimesy;
            if counterxDP = 0 then
                argmin_indexxDN <= counterxDP;
            else
                argmin_indexxDN <= counterxDP - 1;
            end if;
        end if;
        if xtimesy < BETA then
            flipvalidxS <= '1';
        else
            flipvalidxS <= '0';
        end if;
        sum_correlationxDN <= sum_correlationxDP + term3xDP;
        ----------------------------------------------------------------------------------

        case statexDP is

            when IDLE =>

                if startxSI = '1' then

                    --argmin_inverfunc is initially the largest possible value in signed
                    argmin_inverfuncxDN <= (others => '1');
                    argmin_inverfuncxDN(argmin_inverfuncxDN'high) <= '0';

                    counterxDN          <= (others => '0');
                    sum_paritiesxDN     <= (others => '0');
                    sum_correlationxDN  <= (others => '0');
                    term1xDN            <= (others => '0');
                    term2xDN            <= (others => '0');
                    term3xDN            <= (others => '0');
                    parityxDN           <= '1';

                    statexDN <= STEP_ONE;
                
                end if;

            when STEP_ONE =>
                -- common calculation for Objective function, parity check and Inversion function
                for i in 0 to N_BG - 1 loop
                    a1 := BG(to_integer(shift_right(counterxDP,Z_pow)), i); -- BG value at the current position

                    if a1 =- 1 then
                        line_vec1(i) := '0';
                    else
                        a2  := unsigned(a1(Z_bits-2 downto 0));
                        b1  := counterxDP(Z_bits-2 downto 0);
                        line_vec1(i) := xhatxDI(to_integer(a2 + b1) + i*Z);
                    end if;

                end loop;

                vec1_xor := xor line_vec1;
                VecOnexDN(to_integer(counterxDP)) <= vec1_xor;
                parityxDN <= (not vec1_xor) and parityxDP; 


                -- counter logic
                counterxDN <= counterxDP + 1;
                if counterxDP = M - 1 then
                    statexDN   <= STEP_TWO;
                    counterxDN <= (others => '0');

                    -- allows us to exit if parity check is respected
                    if parityxDP = '1' and vec1_xor = '0' then
                        validxDN <= '1';
                        statexDN <= IDLE;
                    end if;
                end if;
            when STEP_TWO =>

                -- Parity term calculation
                if counterxDP < M then
                    if VecOnexDP(to_integer(counterxDP)) = '1' then
                        sum_paritiesxDN <= sum_paritiesxDP - 1;
                    else
                        sum_paritiesxDN <= sum_paritiesxDP + 1;
                    end if;
                end if;

                -- convert y to the t_metric format
                y_resized_shifted := shift_left(resize(yxDI(to_integer(counterxDP)), xtimesy'length), num_metric_fractional_bits - num_fractional_bits);

                if xhatxDI(to_integer(counterxDP)) = '1' then
                    term3xDN <= y_resized_shifted;
                else
                    term3xDN <= - y_resized_shifted;
                end if;

                for i in 0 to M_BG -1 loop
                    a1 := BG(i, to_integer(shift_right(counterxDP,Z_pow))); -- BG value at the current position
                    if a1 /= - 1 then

                        a2 := unsigned(a1(Z_bits-2 downto 0));
                        b1 := counterxDP(Z_bits-2 downto 0);

                        line_vec2(i)  := VecOnexDP(to_integer(b1-a2) + i * Z);
                        line_vec3(i)  := '1';
                    else 
                        line_vec2(i)  := '0';
                        line_vec3(i)  := '0';
                    end if;

                end loop;

                -- sum all bits in line_vec2
                sum_vec2 := (others => '0');
                for i in 0 to M_BG - 1 loop
                    sum_vec2 := sum_vec2 + line_vec2(i);
                end loop;

                -- sum all the bits in (line_vec3 xor line_vec2)
                sum_vec3 := (others => '0');
                for i in 0 to M_BG - 1 loop
                    sum_vec3 := sum_vec3 + (line_vec3(i) xor line_vec2(i));
                end loop;
                
                -- store the values that will be summed in the next pipeline stage
                previndexxDN <= counterxDP;
                term1xDN <= sum_vec2;
                term2xDN <= sum_vec3;

                -- counter logic
                counterxDN <= counterxDP + 1;
                if counterxDP = N - 1 then
                    statexDN <= END_DEC;
                end if;
            
            -- used to delay the valid signal by 1 so that the last pipeline stage can be calculated
            when END_DEC =>
                validxDN <= '1';
                statexDN <= IDLE;

            when others =>
                statexDN <= IDLE;
        end case;

    end process p_main;


    -- valid and parity ok signals
    validxSO  <= validxDP;
    parityxSO <= parityxDP;

    -- the objective function is the sum of the parities and the correlation sums
    valuexDO <= sum_correlationxDP + shift_left(resize(sum_paritiesxDP, sum_correlationxDP'length), num_metric_fractional_bits);

    -- flipping indices
    flipvalidxSO       <= flipvalidxS;
    flipindexxDO       <= previndexxDP;
    SingleFlipindexxDO <= argmin_indexxDP;  -- single bit flip index (valid when validxSO is high)

end architecture rtl;