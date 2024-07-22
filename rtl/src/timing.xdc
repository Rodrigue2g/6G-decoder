create_clock -period 10.000 -name CLKxCI -waveform {0.000 5.000} [get_ports CLKxCI]
set_input_delay -clock [get_clocks CLKxCI] -min -add_delay 0.000 [get_ports {yxDI[*]}]
set_input_delay -clock [get_clocks CLKxCI] -max -add_delay 1.000 [get_ports {yxDI[*]}]
set_input_delay -clock [get_clocks CLKxCI] -min -add_delay 0.000 [get_ports RSTxRI]
set_input_delay -clock [get_clocks CLKxCI] -max -add_delay 1.000 [get_ports RSTxRI]
set_input_delay -clock [get_clocks CLKxCI] -min -add_delay 0.000 [get_ports enablexSI]
set_input_delay -clock [get_clocks CLKxCI] -max -add_delay 1.000 [get_ports enablexSI]
set_output_delay -clock [get_clocks CLKxCI] -min -add_delay 0.000 [get_ports receive_readyxSO]
set_output_delay -clock [get_clocks CLKxCI] -max -add_delay 0.000 [get_ports receive_readyxSO]
set_output_delay -clock [get_clocks CLKxCI] -min -add_delay 0.000 [get_ports valid_outputxSO]
set_output_delay -clock [get_clocks CLKxCI] -max -add_delay 0.000 [get_ports valid_outputxSO]
set_output_delay -clock [get_clocks CLKxCI] -min -add_delay 0.000 [get_ports xhatxDO]
set_output_delay -clock [get_clocks CLKxCI] -max -add_delay 0.000 [get_ports xhatxDO]





