library ieee;
use ieee.std_logic_1164.all;
use ieee.numeric_std.all;
library work;
use work.constants.all;
entity sc_ldpc_decoder_top_level is
    port
    (
        CLKxCI : in std_logic;
        RSTxRI : in std_logic;

        -- input and output (bitstream)
        yxDI    : in signed(FXP_WIDTH - 1 downto 0);
        xhatxDO : out std_logic;

        -- control signals
        enablexSI        : in std_logic;
        receive_readyxSO : out std_logic;
        valid_outputxSO  : out std_logic
    );
end entity sc_ldpc_decoder_top_level;

architecture rtl of sc_ldpc_decoder_top_level is

    

    -- data 
    signal xhatxS                         : std_logic;
    signal y_buffxDP, y_buffxDN           : t_y_array_half;
    signal x_buffxDP, x_buffxDN, x_buffxS : t_xhat_array_half;

    -- counters
    signal y_counterxDP, y_counterxDN : unsigned(N_bits - 1 downto 0);
    signal x_counterxDP, x_counterxDN : unsigned(N_bits - 1 downto 0);

    -- control signals
    signal start_decodingxDP, start_decodingxDN                   : std_logic;
    signal valid_outputxDP, valid_outputxDN                       : std_logic;
    signal receive_readyxDP, receive_readyxDN                     : std_logic;
    signal LDPC_finishedxS                                        : std_logic;
    signal LDPC_readyxS                                           : std_logic;

    -- LDPC decoder 
    component ldpc_decoder_top_level is
        port
        (
            CLKxCI : in std_logic;
            RSTxRI : in std_logic;

            -- input and output (vector) ports
            yxDI    : in t_y_array_half;
            xhatxDO : out t_xhat_array_half;

            -- control signals
            startxSI        : in std_logic;
            input_readyxSO  : out std_logic;
            output_readyxSO : out std_logic

        );
    end component ldpc_decoder_top_level;
begin

    p_dff : process (CLKxCI, RSTxRI)
    begin
        if RSTxRI = '1' then

            -- reset the buffers and counters
            y_counterxDP <= (others => '0');
            y_buffxDP    <= (others => (others => '0'));
            x_counterxDP <= (others => '0');
            x_buffxDP    <= (others => '0');

            -- reset the control signals
            start_decodingxDP <= '0';
            receive_readyxDP  <= '0';
            valid_outputxDP   <= '0';
        elsif rising_edge(CLKxCI) then

            -- update the buffers
            y_counterxDP <= y_counterxDN;
            y_buffxDP    <= y_buffxDN;
            x_counterxDP <= x_counterxDN;
            x_buffxDP    <= x_buffxDN;

            -- update the control signals
            start_decodingxDP <= start_decodingxDN;
            receive_readyxDP  <= receive_readyxDN;
            valid_outputxDP   <= valid_outputxDN;

        end if;
    end process p_dff;

    -- controls the signals for the input and LDPC decoder
    p_input_control : process (all)
    begin

        -- default values
        y_counterxDN      <= y_counterxDP;
        start_decodingxDN <= '0';
        receive_readyxDN  <= '1';

        if enablexSI = '1' then

            y_counterxDN <= y_counterxDP + 1;

            -- if the buffer is full and the LDPC is not running we can start the decoding process
            if y_counterxDP = to_unsigned(N/2 - 1, N_bits) then

                if LDPC_readyxS = '1' then
                    y_counterxDN      <= (others => '0');
                    start_decodingxDN <= '1';
                    receive_readyxDN  <= '1';

                else
                    y_counterxDN      <= to_unsigned(N/2 - 1, y_counterxDN'length);
                    start_decodingxDN <= '0';
                    receive_readyxDN  <= '0';

                end if;
            end if;
        else
            y_counterxDN <= (others => '0');
        end if;
    end process p_input_control;


    -- -- loads the input buffer with the values
    p_input_buffer : process (all)
    begin
        if receive_readyxDP = '1' then
            -- fill the input buffer with the value and increment the buffer position counter
            for i in 0 to (N/2 - 2) loop
                y_buffxDN(i) <= y_buffxDP(i + 1);
            end loop;
            y_buffxDN(N/2 - 1) <= yxDI;
        else
            y_buffxDN <= y_buffxDP;
        end if;
    end process p_input_buffer;


    p_output : process (all)
    begin

        -- default values
        x_counterxDN    <= x_counterxDP;
        valid_outputxDN <= valid_outputxDP;
        xhatxS <= '0';

        x_buffxDN <= x_buffxDP;
        if valid_outputxDP = '1' then
            x_counterxDN <= x_counterxDP + 1;

            if x_counterxDP = to_unsigned(N/2 - 1, N_bits) then
                valid_outputxDN <= '0';
            end if;

            -- shift and push the value out
            for i in 0 to (N/2 - 2) loop
                x_buffxDN(i) <= x_buffxDP(i+1);
            end loop;
            xhatxS  <= x_buffxDP(0);
        end if;

        -- if the LDPC has finished we can start outputing the values
        if LDPC_finishedxS = '1' then
            valid_outputxDN <= '1';
            x_counterxDN    <= (others => '0');
            x_buffxDN       <= x_buffxS;
        end if;

    end process p_output;

    ldpc_decoder : ldpc_decoder_top_level
    port map
    (
        CLKxCI => CLKxCI,
        RSTxRI => RSTxRI,

        yxDI    => y_buffxDP,
        xhatxDO => x_buffxS,

        startxSI        => start_decodingxDP,
        output_readyxSO => LDPC_finishedxS,
        input_readyxSO  => LDPC_readyxS
    );

    receive_readyxSO <= receive_readyxDP;
    valid_outputxSO  <= valid_outputxDP;
    xhatxDO          <= xhatxS;

end architecture rtl;