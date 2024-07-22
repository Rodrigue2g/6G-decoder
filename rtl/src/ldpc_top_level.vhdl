library ieee;
use ieee.std_logic_1164.all;
use ieee.numeric_std.all;


library work;
use work.constants.all;


entity ldpc_decoder_top_level is 
    port (
        CLKxCI : in std_logic;
        RSTxRI : in std_logic;

        -- input and output (vector) ports
        yxDI    : in  t_y_array_half;
        xhatxDO : out t_xhat_array_half;

        -- control signals
        startxSI            : in std_logic;
        input_readyxSO      : out std_logic;
        output_readyxSO     : out std_logic
    );
end entity ldpc_decoder_top_level;

architecture rtl of ldpc_decoder_top_level is

    type t_state is (IDLE, STEP_ONE, STEP_TWO, STEP_THREE, STEP_FOUR, STEP_FIVE, STEP_WAIT);

    -- ! might change the logic
    signal first_runxDP, first_runxDN : std_logic;

    -- signals
    signal startMetricxDP, startMetricxDN : std_logic;
    signal modexDN, modexDP : std_logic;
    signal MetricValidxS, ParCheckParityxS : std_logic;
    signal ObjFuncValue_lastxDN, ObjFuncValue_lastxDP : t_metric;
    signal statexDN, statexDP : t_state;
    signal ObjFuncValuexS: t_metric;
    signal flipvalidxS : std_logic;
    signal flipindexxS : unsigned(N_bits-1 downto 0);
    signal SingleFlipindexxS : unsigned(N_bits-1 downto 0);
    signal ArmReadyxDN, ArmReadyxDP : std_logic_vector(2-1 downto 0);
    signal output_counterxDP, output_counterxDN : unsigned(N_bits-1 downto 0);

    -- buffers and counters
    signal y_buffxDP, y_buffxDN : t_y_array;
    signal x_buffxDP, x_buffxDN : t_xhat_array;
    signal bitmaskxDP, bitmaskxDN : t_xhat_array;
    signal iterxDP, iterxDN : unsigned(MAX_ITERATIONS_bits-1 downto 0);

    component MetricCalc is
        port (
            CLKxCI : in std_logic;
            RSTxRI : in std_logic;
    
            xhatxDI : in t_xhat_array;
            yxDI    : in t_y_array;
    
            startxSI: in std_logic;
    
            validxSO    : out std_logic;
            parityxSO   : out std_logic;
            valuexDO    : out t_metric;

            flipvalidxSO: out std_logic;
            flipindexxDO: out unsigned(N_bits-1 downto 0);

            SingleFlipindexxDO: out unsigned(N_bits-1 downto 0)
        );
    end component MetricCalc;

    begin  

        p_dff : process(CLKxCI, RSTxRI)
        begin
            if (RSTxRI = '1') then

                -- buffers and counters
                y_buffxDP <= (others => (others => '0'));
                x_buffxDP <= (others => '0');
                iterxDP <= (others => '0');
                bitmaskxDP <= (others => '0');

                -- signals
                startMetricxDP <= '0';
                modexDP <= '0';
                ObjFuncValue_lastxDP <= (others => '0');
                statexDP <= IDLE;
                first_runxDP <= '1';
                ArmReadyxDP <= (others => '0');
                output_counterxDP <= (others => '0');

            elsif rising_edge(CLKxCI) then

                -- buffers and counters
                y_buffxDP <= y_buffxDN;
                x_buffxDP <= x_buffxDN;
                iterxDP <= iterxDN;
                bitmaskxDP <= bitmaskxDN;

                -- signals
                startMetricxDP <= startMetricxDN;
                modexDP <= modexDN;
                ObjFuncValue_lastxDP <= ObjFuncValue_lastxDN;
                statexDP <= statexDN;
                first_runxDP <= first_runxDN;
                ArmReadyxDP <= ArmReadyxDN;
                output_counterxDP <= output_counterxDN;

            end if;
        end process p_dff;


        p_dec : process(all)
        begin

            -- default values
            y_buffxDN <= y_buffxDP;
            x_buffxDN <= x_buffxDP;
            startMetricxDN <= '0';
            first_runxDN <= first_runxDP;
            modexDN <= modexDP;
            ObjFuncValue_lastxDN <= ObjFuncValue_lastxDP;
            statexDN <= statexDP;
            iterxDN <= iterxDP;
            bitmaskxDN <= bitmaskxDP;
            ArmReadyxDN <= ArmReadyxDP;

            -- this counter makes sure that the output isn't pushed out before the previous output has been outputed
            output_counterxDN <= output_counterxDP+1;
            if output_counterxDP = to_unsigned(N/2, output_counterxDP'length) then
                ArmReadyxDN(0) <= '1';
            end if;

            case statexDP is
                when IDLE =>

                    if startxSI = '1' then

                        -- set to zero the first run flag
                        first_runxDN <= '0';

                        -- start decoding
                        statexDN <= IDLE      when first_runxDP = '1' else STEP_ONE;
                        startMetricxDN <= '0' when first_runxDP = '1' else '1';
                        
                        -- (1,0) -> (-1,1) -> (0,1)
                        -- shift and fill the buffer
                        for i in 0 to N/2-1 loop
                            x_buffxDN(i)      <= x_buffxDP(i+N/2);
                            y_buffxDN(i)      <= y_buffxDP(i+N/2);
                        end loop;
                        for i in N/2 to N-1 loop
                            x_buffxDN(i)      <= '1' when yxDI(i-N/2)>=0 else '0';
                            y_buffxDN(i)      <= yxDI(i-N/2);
                        end loop;       
                            
                        -- start in multi-bit mode
                        modexDN <= '0';
                        ObjFuncValue_lastxDN <= (others => '0');
                        iterxDN <= (others => '0');
                        bitmaskxDN <= (others => '0');
                        output_counterxDN <= (others => '0');
                        ArmReadyxDN <= (others => '0');
                    end if;
                when STEP_ONE =>

                    -- flip the bits if we have a valid flip index
                    if modexDP = '0' then
                        if flipvalidxS = '1' then
                            bitmaskxDN(to_integer(flipindexxS)) <= '1';
                        end if;
                    elsif MetricValidxS='1' then
                        bitmaskxDN(to_integer(SingleFlipindexxS)) <= '1';
                    end if;   

                    -- if metric finished
                    if MetricValidxS = '1' then

                        -- if parity check is ok, we have converged
                        if ParCheckParityxS = '1' then

                            ArmReadyxDN(1) <= '1';
                            statexDN <= STEP_WAIT;

                        else
                            -- else we remember the current obj value
                            ObjFuncValue_lastxDN <= ObjFuncValuexS;
                            statexDN <= STEP_TWO;
                        end if;
                    end if; 

                when STEP_TWO =>
                    
                    -- flip all the bits in x from the bitmask
                    statexDN <= STEP_THREE;
                    startMetricxDN <= '1';
                    for i in 0 to N-1 loop
                        x_buffxDN(i) <= x_buffxDP(i) xor bitmaskxDP(i);
                    end loop;
                    

                when STEP_THREE=>
                    
                    if MetricValidxS = '1' then

                        -- increment the iteration counter
                        iterxDN <= iterxDP + 1;

                        -- if the parity check is ok, we have converged
                        if ParCheckParityxS = '1' or iterxDP = to_unsigned(MAX_ITERATIONS-1, iterxDP 'length) then

                            ArmReadyxDN(1) <= '1';
                            statexDN <= STEP_WAIT;

                        else
                            statexDN <= STEP_FOUR;

                            -- if the last obj value is greater than the current obj value
                            if ObjFuncValue_lastxDP >= ObjFuncValuexS then
                                
                                -- if we are in multi-bit mode we switch to single bit mode
                                if modexDP = '0' then
                                    modexDN <= '1'; 
                                
                                -- if we are in single bit mode, we have converged
                                else 
                                    -- exit
                                    statexDN <= STEP_FIVE;
                                end if;
                            else

                                -- reset the bit mask for the next iteration
                                bitmaskxDN <= (others => '0');

                            end if;        
                        end if;
                    end if;
                
                when STEP_FOUR =>
                    statexDN <= STEP_ONE;
                    startMetricxDN <= '1';

                    for i in 0 to N-1 loop
                        x_buffxDN(i) <= x_buffxDP(i) xor bitmaskxDP(i);
                    end loop;

                    bitmaskxDN <= (others => '0');

                when STEP_FIVE =>
                    
                    for i in 0 to N-1 loop
                        x_buffxDN(i) <= x_buffxDP(i) xor bitmaskxDP(i);
                    end loop;

                    ArmReadyxDN(1) <= '1';
                    statexDN <= STEP_WAIT;
                    
                -- this step is to make sure that the previous output has been fully outputed
                when STEP_WAIT =>
                    if ArmReadyxDP = "11" then
                        statexDN <= IDLE;
                    end if;

                when others =>
                    statexDN <= IDLE;
            end case;

            -- reset the output ready signal
            if ArmReadyxDP = "11" then
                ArmReadyxDN <= (others => '0');
            end if;
           
      
        end process p_dec;



        -- instantiate the parity check module
        metric_calculation : MetricCalc
            port map (
                CLKxCI => CLKxCI,
                RSTxRI => RSTxRI,
                yxDI => y_buffxDP,
                xhatxDI => x_buffxDP,
                startxSI => startMetricxDP,
                validxSO => MetricValidxS,
                parityxSO => ParCheckParityxS,
                valuexDO => ObjFuncValuexS,
                flipvalidxSO => flipvalidxS,
                flipindexxDO => flipindexxS,
                SingleFlipindexxDO => SingleFlipindexxS
            );


        -- output signals
        output_readyxSO <= '1' when ArmReadyxDP="11" else '0';  -- if the output is valid
        input_readyxSO  <= '1' when statexDP = IDLE else '0';   -- if we are ready to accept new data
            
        process(all)
        begin
            for i in 0 to N/2-1 loop
                xhatxDO(i) <= x_buffxDP(i);                     -- output the first half of the codeword
            end loop;
        end process;

  


end architecture rtl;
    