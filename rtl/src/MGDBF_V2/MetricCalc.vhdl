library ieee;
use ieee.std_logic_1164.all;
use ieee.numeric_std.all;

library work;
use work.constants.all;

-- ! testbench
--use ieee.std_logic_textio.all;
-- Takes 128 cycles max
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
    type t_state is (IDLE, STEP_ONE, STEP_TWO);

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
    signal flipvalidxDN, flipvalidxDP               : std_logic;
    signal flipindexxDN, flipindexxDP               : unsigned(N_bits - 1 downto 0);
    signal argmin_inverfuncxDN, argmin_inverfuncxDP : t_metric;
    signal argmin_indexxDN, argmin_indexxDP         : unsigned(N_bits - 1 downto 0);
    signal statexDN, statexDP                       : t_state;

begin

    BG <= (
        (to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(13, Z_bits), to_signed(-1, Z_bits), to_signed(14, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(2, Z_bits), to_signed(5, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(9, Z_bits), to_signed(3, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(4, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(3, Z_bits), to_signed(-1, Z_bits), to_signed(13, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(13, Z_bits), to_signed(-1, Z_bits), to_signed(11, Z_bits), to_signed(7, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(9, Z_bits), to_signed(14, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(9, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(6, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(14, Z_bits), to_signed(2, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(0, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits)),
        (to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(12, Z_bits), to_signed(-1, Z_bits), to_signed(6, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(1, Z_bits), to_signed(8, Z_bits), to_signed(-1, Z_bits), to_signed(14, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(13, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(10, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(4, Z_bits), to_signed(-1, Z_bits), to_signed(8, Z_bits), to_signed(4, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(5, Z_bits), to_signed(-1, Z_bits), to_signed(2, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(11, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(6, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(10, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(12, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(4, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(12, Z_bits), to_signed(-1, Z_bits), to_signed(0, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits)),
        (to_signed(14, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(9, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(0, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(2, Z_bits), to_signed(-1, Z_bits), to_signed(10, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(0, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(15, Z_bits), to_signed(13, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(8, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(3, Z_bits), to_signed(5, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(14, Z_bits), to_signed(9, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(2, Z_bits), to_signed(-1, Z_bits), to_signed(0, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(8, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(5, Z_bits), to_signed(0, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(3, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(0, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits)),
        (to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(3, Z_bits), to_signed(-1, Z_bits), to_signed(10, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(11, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(12, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(9, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(10, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(0, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(3, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(11, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(1, Z_bits), to_signed(1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(5, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(11, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(2, Z_bits), to_signed(15, Z_bits), to_signed(-1, Z_bits), to_signed(8, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(12, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(2, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(3, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(0, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits)),
        (to_signed(5, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(3, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(11, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(7, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(5, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(5, Z_bits), to_signed(-1, Z_bits), to_signed(5, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(9, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(2, Z_bits), to_signed(1, Z_bits), to_signed(1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(13, Z_bits), to_signed(13, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(15, Z_bits), to_signed(9, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(6, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(3, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(1, Z_bits), to_signed(13, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(0, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits)),
        (to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(5, Z_bits), to_signed(1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(6, Z_bits), to_signed(-1, Z_bits), to_signed(3, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(14, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(13, Z_bits), to_signed(-1, Z_bits), to_signed(15, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(14, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(2, Z_bits), to_signed(13, Z_bits), to_signed(12, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(5, Z_bits), to_signed(1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(7, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(11, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(0, Z_bits), to_signed(6, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(13, Z_bits), to_signed(-1, Z_bits), to_signed(7, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(0, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits)),
        (to_signed(-1, Z_bits), to_signed(11, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(6, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(14, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(12, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(13, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(10, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(10, Z_bits), to_signed(-1, Z_bits), to_signed(4, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(11, Z_bits), to_signed(9, Z_bits), to_signed(-1, Z_bits), to_signed(13, Z_bits), to_signed(-1, Z_bits), to_signed(11, Z_bits), to_signed(-1, Z_bits), to_signed(1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(7, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(1, Z_bits), to_signed(1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(4, Z_bits), to_signed(-1, Z_bits), to_signed(5, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(0, Z_bits), to_signed(-1, Z_bits)),
        (to_signed(-1, Z_bits), to_signed(7, Z_bits), to_signed(10, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(7, Z_bits), to_signed(2, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(0, Z_bits), to_signed(11, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(15, Z_bits), to_signed(-1, Z_bits), to_signed(11, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(9, Z_bits), to_signed(8, Z_bits), to_signed(11, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(0, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(13, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(10, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(14, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(9, Z_bits), to_signed(3, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(12, Z_bits), to_signed(-1, Z_bits), to_signed(11, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(-1, Z_bits), to_signed(0, Z_bits))
        );

    p_dff : process (CLKxCI, RSTxRI)
    begin
        if RSTxRI = '1' then
            validxDP            <= '0';
            parityxDP           <= '0';
            flipvalidxDP        <= '0';
            flipindexxDP        <= (others => '0');
            counterxDP          <= (others => '0');
            sum_paritiesxDP     <= (others => '0');
            sum_correlationxDP  <= (others => '0');
            VecOnexDP           <= (others => '0');
            argmin_inverfuncxDP <= (others => '0');
            argmin_indexxDP     <= (others => '0');
            statexDP            <= IDLE;

        elsif rising_edge(CLKxCI) then
            validxDP            <= validxDN;
            parityxDP           <= parityxDN;
            counterxDP          <= counterxDN;
            flipvalidxDP        <= flipvalidxDN;
            flipindexxDP        <= flipindexxDN;
            sum_paritiesxDP     <= sum_paritiesxDN;
            sum_correlationxDP  <= sum_correlationxDN;
            VecOnexDP           <= VecOnexDN;
            argmin_inverfuncxDP <= argmin_inverfuncxDN;
            argmin_indexxDP     <= argmin_indexxDN;
            statexDP            <= statexDN;
        end if;
    end process p_dff;
    p_main : process (all)
        variable line_vec1 : std_logic_vector(0 to N_BG - 1);
        variable xtimesy   : t_metric;
        variable a, b      : integer;

        variable extended_VecOnexDP : t_metric;

    begin

        -- default values
        validxDN            <= '0';
        parityxDN           <= parityxDP;
        counterxDN          <= counterxDP;
        flipvalidxDN        <= '0';
        flipindexxDN        <= (others => '0');
        sum_paritiesxDN     <= sum_paritiesxDP;
        sum_correlationxDN  <= sum_correlationxDP;
        VecOnexDN           <= VecOnexDP;
        argmin_inverfuncxDN <= argmin_inverfuncxDP;
        argmin_indexxDN     <= argmin_indexxDP;
        statexDN            <= statexDP;

        case statexDP is

            when IDLE =>

                if startxSI = '1' then
                    counterxDN <= (others => '0');

                    sum_paritiesxDN    <= (others => '0');
                    sum_correlationxDN <= (others => '0');
                    parityxDN          <= '1';

                    VecOnexDN                                     <= (others => '0');
                    argmin_inverfuncxDN                           <= (others => '1');
                    argmin_inverfuncxDN(argmin_inverfuncxDN'high) <= '0'; -- set the sign bit to 0
                    argmin_indexxDN                               <= (others => '0');

                    statexDN <= STEP_ONE;
                end if;

            when STEP_ONE =>

                -- common calculation for Objective function, parity check and Inversion function
                for i in 0 to N_BG - 1 loop
                    a := to_integer(BG(to_integer(counterxDP/Z), i)); -- BG value at the current position

                    if a =- 1 then
                        line_vec1(i) := '0';
                    else
                        b := to_integer(counterxDP) mod Z;

                        line_vec1(i) := xhatxDI((a + b) mod Z + i * Z);
                    end if;

                end loop;

                VecOnexDN(to_integer(counterxDP)) <= xor line_vec1; --! Check xor implementation (long critical path ?)
                -- counter logic
                counterxDN <= counterxDP + 1;
                if counterxDP = M - 1 then
                    statexDN   <= STEP_TWO;
                    counterxDN <= (others => '0');
                end if;
            when STEP_TWO =>

                -- Objective function + parity check
                if counterxDP < M then
                    if VecOnexDP(to_integer(counterxDP)) = '1' then
                        parityxDN       <= '0';
                        sum_paritiesxDN <= sum_paritiesxDP - 1;
                    else
                        sum_paritiesxDN <= sum_paritiesxDP + 1;
                    end if;
                end if;

                if xhatxDI(to_integer(counterxDP)) = '1' then
                    xtimesy := shift_left(resize(yxDI(to_integer(counterxDP)), xtimesy'length), num_metric_fractional_bits - num_fractional_bits);
                else
                    xtimesy := - shift_left(resize(yxDI(to_integer(counterxDP)), xtimesy'length), num_metric_fractional_bits - num_fractional_bits);
                end if;
                sum_correlationxDN <= sum_correlationxDP + xtimesy;

                -- ! probably long critical path ?
                -- Inversion function, xtimesy = xtimesy is the result of the inversion function for the current bit
                for i in 0 to M_BG - 1 loop
                    a := to_integer(BG(i, to_integer(counterxDP/Z))); -- BG value at the current position

                    if a /= - 1 then

                        b := to_integer(counterxDP) mod Z;

                        extended_VecOnexDP                             := (others => '0');
                        extended_VecOnexDP(num_metric_fractional_bits) := '1';

                        if VecOnexDP((Z - a + b) mod Z + i * Z) = '1' then
                            xtimesy := xtimesy - extended_VecOnexDP;

                        else
                            xtimesy := xtimesy + extended_VecOnexDP;
                        end if;
                    end if;

                end loop;

                -- Single bit flip logic
                if argmin_inverfuncxDP > xtimesy then
                    argmin_inverfuncxDN <= xtimesy;
                    argmin_indexxDN     <= counterxDP;
                end if;

                -- Multi-bit flip logic
                if xtimesy < BETA then
                    flipvalidxDN <= '1';
                    flipindexxDN <= counterxDP;
                end if;

                -- counter logic
                counterxDN <= counterxDP + 1;
                if counterxDP = N - 1 then
                    statexDN <= IDLE;
                    validxDN <= '1';
                end if;

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
    flipvalidxSO       <= flipvalidxDP; -- multibit flip index valid
    flipindexxDO       <= flipindexxDP; -- multibit flip index
    SingleFlipindexxDO <= argmin_indexxDP; -- single bit flip index (valid when validxSO is high)
end architecture rtl;