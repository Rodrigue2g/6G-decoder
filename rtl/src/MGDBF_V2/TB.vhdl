library ieee;
use ieee.std_logic_1164.all;
use ieee.numeric_std.all;
use ieee.math_real.all;

-- use sfixed
use ieee.fixed_pkg.all;
use std.textio.all;

library work;
use work.constants.all;

entity TB is
end entity TB;
architecture rtl of TB is
    constant CLK_HIGH : time := 4 ns;
    constant CLK_LOW  : time := 4 ns;
    constant CLK_PER  : time := CLK_LOW + CLK_HIGH;
    constant CLK_STIM : time := 1 ns;

    signal clk : std_logic := '0';
    signal rst : std_logic := '0';

    -- change the path for the input file
    file input_file  : text open read_mode  is "/home/mike/Documents/EPFL/BA6/PROJET/24s-michel-rodrigue-bitflippingdecoder6g/VHDL/MGDBF_V2/y_out.txt";
    file output_file : text open write_mode is "/home/mike/Documents/EPFL/BA6/PROJET/24s-michel-rodrigue-bitflippingdecoder6g/VHDL/MGDBF_V2/xhat_out.txt";
    file verif_file  : text open read_mode  is "/home/mike/Documents/EPFL/BA6/PROJET/24s-michel-rodrigue-bitflippingdecoder6g/VHDL/MGDBF_V2/xhat_out_MATLAB.txt";

    signal y             : signed(FXP_WIDTH - 1 downto 0) := (others => '0');
    signal xhat          : std_logic                      := '0';
    signal enable        : std_logic                      := '0';
    signal receive_ready : std_logic;
    signal valid_output  : std_logic;

    signal error : std_logic := '0';

    -- DUT component
    component sc_ldpc_decoder_top_level is
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
    end component sc_ldpc_decoder_top_level;

begin

    -- Clock generation
    p_clock : process is
    begin
        clk <= '0';
        wait for CLK_LOW;
        clk <= '1';
        wait for CLK_HIGH;
    end process p_clock;

    -- Reset generation
    p_reset : process is
    begin
        rst <= '1';
        wait until clk'event and clk = '1'; -- wait for rising edge
        wait for (1 * CLK_PER + CLK_STIM);
        rst <= '0';
        wait;
    end process p_reset;

    -- y process (reads from file and sends to DUT)
    y_generation : process is
        variable line        : line;
        variable real_number : real;
        variable rounded_num : real;
    begin

        -- delay start
        wait until rst = '1';
        wait for 2 * CLK_PER;
        wait until clk'event and clk = '1';

        -- loop while not end of file
        while not endfile(input_file) loop

            if receive_ready = '1' then
                -- read a new value from the file and converts the real number to a fixed point number
                readline(input_file, line);
                read(line, real_number);

                rounded_num := real_number;
                y      <= signed(to_sfixed(rounded_num * scale_factor, y'high, y'low));
                enable <= '1'; -- enable the sc-ldpc at first data input
            end if;

            wait for CLK_PER;

        end loop;
        -- close the file and stop the simulation
        file_close(input_file);

        y <= (others => '0');

        wait;

    end process y_generation;

    -- xhat process (writes to file)
    xhat_generation : process is
        variable line_write    : line;
        variable line_read     : line;
        variable counts        : integer := 0;
        variable xhat_out_real : real    := 0.0;

    begin

        -- delay start
        wait until rst = '1';
        wait for 2 * CLK_PER;
        wait until clk'event and clk = '1';

        -- loop while not end of file (32000 - last codeword)
        while counts < 32000 - 640 loop

            error <= '0';
            if valid_output = '1' then
                -- write the value to the file
                write(line_write, xhat);
                writeline(output_file, line_write);

                readline(verif_file, line_read);
                read(line_read, xhat_out_real);

                if xhat_out_real < 0.5 and xhat = '1' then
                    error <= '0';
                elsif xhat_out_real >= 0.5 and xhat = '0' then
                    error <= '0';
                else
                    error <= '1';
                end if;

                counts := counts + 1;
            end if;

            wait for CLK_PER;

        end loop;

        -- close the file and stop the simulation
        file_close(output_file);
        file_close(verif_file);

        -- stop the simulation 
        assert false report "Whole message decoded" severity failure;

    end process xhat_generation;
    
    -- DUT instantiation
    DUT_SC_LDPC_top : sc_ldpc_decoder_top_level
    port map
    (
        CLKxCI  => clk,
        RSTxRI  => rst,
        yxDI    => y,
        xhatxDO => xhat,

        enablexSI        => enable,
        receive_readyxSO => receive_ready,
        valid_outputxSO  => valid_output
    );
end architecture rtl;