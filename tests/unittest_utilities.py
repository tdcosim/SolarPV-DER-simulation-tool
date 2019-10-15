"""Common utility functions used in unit tests."""



def show_DER_status(pvder_object):
        """Show DER states."""     
        
        pvder_object.show_PV_DER_states(quantity='power')
        pvder_object.show_PV_DER_states(quantity='current')
        pvder_object.show_PV_DER_states(quantity='voltage')                
        pvder_object.show_PV_DER_states(quantity='duty cycle')
        pvder_object.show_RT_settings(settings_type='LVRT')
        pvder_object.show_RT_settings(settings_type='HVRT')   
    
def plot_DER_trajectories(results_object):
        """PLot DER trajectory."""
        
        results_object.PER_UNIT = False
        results_object.PLOT_TITLE = True
        results_object.font_size = 18
        
        results_object.plot_DER_simulation(plot_type='active_power_Ppv_Pac_PCC')#
        results_object.plot_DER_simulation(plot_type='reactive_power')
        results_object.plot_DER_simulation(plot_type='current')
        results_object.plot_DER_simulation(plot_type='voltage_Vdc')
        results_object.plot_DER_simulation(plot_type='voltage_LV')
        results_object.plot_DER_simulation(plot_type='voltage_Vpcclv_all_phases')
        results_object.plot_DER_simulation(plot_type='duty_cycle')  