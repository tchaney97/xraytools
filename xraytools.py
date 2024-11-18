import numpy as np
from prettytable import PrettyTable
import matplotlib.pyplot as plt
from matplotlib.pyplot import subplots
import xraydb

def calc_sdd(length_y, length_x, energy, max_q):
    """
    Calculate the sample-to-detector distance (SDD) required to achieve a given maximum scattering vector (Q) at the detector edge or corner.

    Args:
        length_y (float): The vertical length of the detector in mm.
        length_x (float): The horizontal length of the detector in mm.
        energy (float): The energy of the X-ray beam in eV.
        max_q (float): The maximum scattering vector in Å^-1.

    Returns:
        None: Prints a table showing the SDD for different beam positions (center, bottom, side, corner) on the detector.
    """
    wavelength = 12400/energy # in Å
    max_theta = 2*(np.arcsin(max_q*wavelength/(4*np.pi)))
    
    # beam position = center
    min_length = min(length_x/2, length_y/2)
    edge_d = min_length
    corner_d = ((length_y/2)**2 + (length_x/2)**2)**0.5
    center_sdd_e = edge_d/np.tan(max_theta)
    center_sdd_c = sdd_corner = corner_d/np.tan(max_theta)
    # beam position = bottom
    min_length = min(length_x/2, length_y)
    edge_d = min_length
    corner_d = ((length_y)**2 + (length_x/2)**2)**0.5
    bottom_sdd_e = edge_d/np.tan(max_theta)
    bottom_sdd_c = sdd_corner = corner_d/np.tan(max_theta)
    #beam position = side
    min_length = min(length_x, length_y/2)
    edge_d = min_length
    corner_d = ((length_y/2)**2 + (length_x)**2)**0.5
    side_sdd_e = edge_d/np.tan(max_theta)
    side_sdd_c = sdd_corner = corner_d/np.tan(max_theta)
    # beam position = corner
    min_length = min(length_x, length_y)
    edge_d = min_length
    corner_d = ((length_y)**2 + (length_x)**2)**0.5
    corner_sdd_e = edge_d/np.tan(max_theta)
    corner_sdd_c = sdd_corner = corner_d/np.tan(max_theta)
    
    # plot table
    table = PrettyTable()
    table.add_column('', ['SDD (max ring)', 'SDD (corner)'])
    table.add_column('Center Beam', [center_sdd_e, center_sdd_c])
    table.add_column('Bottom Beam', [bottom_sdd_e, bottom_sdd_c])
    table.add_column('Side Beam', [side_sdd_e, side_sdd_c])
    table.add_column('Corner Beam', [corner_sdd_e, corner_sdd_c])
    table.float_format = '.0'
    print(f'Table showing SDDs (in mm) to put Q={max_q:.2f} Å-1 at the edge/corner of the detector')
    print(table)

def calc_qrange(length_y, length_x, bs_radius, energy, sdd):
    """
    Calculate the range of scattering vectors (Q) for a given sample-to-detector distance (SDD).

    Args:
        length_y (float): The vertical length of the detector in mm.
        length_x (float): The horizontal length of the detector in mm.
        bs_radius (float): The radius of the beamstop in mm.
        energy (float): The energy of the X-ray beam in eV.
        sdd (float): The sample-to-detector distance in mm.

    Returns:
        None: Prints a table showing the maximum and minimum Q values for different beam positions (center, bottom, side, corner) on the detector.
    """
    wavelength = 12400/energy # in Å
    
    # beam position = center
    min_length = min(length_x/2, length_y/2)
    edge_d = min_length
    corner_d = ((length_y/2)**2 + (length_x/2)**2)**0.5
    center_q_e = 4*np.pi*np.sin(np.arctan(edge_d/sdd)/2)/wavelength
    center_q_c = 4*np.pi*np.sin(np.arctan(corner_d/sdd)/2)/wavelength
    
    # beam position = bottom
    min_length = min(length_x/2, length_y)
    edge_d = min_length
    corner_d = ((length_y)**2 + (length_x/2)**2)**0.5
    bottom_q_e = 4*np.pi*np.sin(np.arctan(edge_d/sdd)/2)/wavelength
    bottom_q_c = 4*np.pi*np.sin(np.arctan(corner_d/sdd)/2)/wavelength
    
    #beam position = side
    min_length = min(length_x, length_y/2)
    edge_d = min_length
    corner_d = ((length_y/2)**2 + (length_x)**2)**0.5
    side_q_e = 4*np.pi*np.sin(np.arctan(edge_d/sdd)/2)/wavelength
    side_q_c = 4*np.pi*np.sin(np.arctan(corner_d/sdd)/2)/wavelength
    
    # beam position = corner
    min_length = min(length_x, length_y)
    edge_d = min_length
    corner_d = ((length_y)**2 + (length_x)**2)**0.5
    corner_q_e = 4*np.pi*np.sin(np.arctan(edge_d/sdd)/2)/wavelength
    corner_q_c = 4*np.pi*np.sin(np.arctan(corner_d/sdd)/2)/wavelength
    
    # minimum q
    min_q = 4*np.pi*np.sin(np.arctan(bs_radius/sdd)/2)/wavelength
    
    # Convert min_q to scientific notation
    min_q_sci = f"{min_q:.1e}"
    
    # plot table
    table = PrettyTable()
    table.add_column('', ['Qmax (max ring)', 'Qmax (corner)', 'Qmin'])
    table.add_column('Center Beam', [f"{center_q_e:.2f}", f"{center_q_c:.2f}", min_q_sci])
    table.add_column('Bottom Beam', [f"{bottom_q_e:.2f}", f"{bottom_q_c:.2f}", min_q_sci])
    table.add_column('Side Beam', [f"{side_q_e:.2f}", f"{side_q_c:.2f}", min_q_sci])
    table.add_column('Corner Beam', [f"{corner_q_e:.2f}", f"{corner_q_c:.2f}", min_q_sci])
    print(f'Table showing extreme Q values (Å-1) given an SDD of {sdd:.0f} mm')
    print(table)

def calc_critical_angle(energy, stoichiometry, density, verbose=True):
    """
    Calculate the critical angle for total external reflection for a given material at a specified energy.

    Args:
        energy (float): The energy of the X-ray beam in eV.
        stoichiometry (str): The chemical formula of the material.
        density (float): The density of the material in g/cm³.
        verbose (bool, optional): If True, prints the critical angle. Default is True.

    Returns:
        float: The critical angle in degrees.
    """
    delta = xraydb.xray_delta_beta(stoichiometry, density, energy)[0]
    critical_angle = np.rad2deg(np.sqrt(2*delta))
    if verbose:
        print(f'Critical angle for {stoichiometry} at {energy}eV is {critical_angle:0.2f}°')
    
    return critical_angle

def calc_xefi(energy, 
              film_stoichiometry, 
              film_density, 
              sub_stoichiometry, 
              sub_density,
              sampthick, 
              incidentangs,
              plot=False):
    """
    Calculate the X-ray Electric Field Intensity (XEFI) throughout a thin film's depth. Assumes one material \
    and zero surface roughness. 

    Args:
        energy (float): The energy of the X-ray beam in eV.
        film_stoichiometry (str): The chemical formula of the film material.
        film_density (float): The density of the film material in g/cm³.
        sub_stoichiometry (str): The chemical formula of the substrate material.
        sub_density (float): The density of the substrate material in g/cm³.
        sampthick (float): Thickness of thin film in nm 
        incidentang (float): Angle of incidence in degrees
        incidentang_extent) (float): +/- extent around incidentang for range. (i.e. 0.9 ± 2, 0.7 to 0.11 range)

    Returns:
        tuple of numpy arrays: (aois, depth, xefi)
            aois  : dimension/coordinate 1D array - angles of incidence
            depth : dimension/coordinate 1D array - depth
            xefi  : electric field data 2D array
        
    """
    
    delta1, beta1, attlen1 = xraydb.xray_delta_beta(film_stoichiometry, film_density, energy)
    delta2, beta2, attlen2 = xraydb.xray_delta_beta(sub_stoichiometry, sub_density, energy)

    beamdiverge = 0.01  # beam divergence in degrees
    
    # Constants
    h = 6.626e-34      # Planck's constant (J*s)
    c = 3e8            # Speed of light (m/s)
    eV_to_J = 1.602e-19  # Conversion factor, eV to J
    lambda_angstrom = (h * c) / (energy * eV_to_J) * 1e10  # Wavelength in Angstroms
    sampthick_angstrom = sampthick * 10  # Convert thickness from nm to Angstroms
    
    angresamp = incidentangs * np.pi / 180
    
    # Define wavevectors
    k0 = 2 * np.pi / lambda_angstrom
    kz0 = k0 * np.sin(angresamp)
    kc1 = np.sqrt(2 * delta1) * k0
    kc2 = np.sqrt(2 * delta2) * k0
    
    # Wavevectors in film (kz1) and substrate (kz2)
    p1 = 1/np.sqrt(2) * np.sqrt(np.sqrt((kz0**2-kc1**2)**2 + 4 * beta1**2 * k0**4) - kc1**2 + kz0**2)
    q1 = 1/np.sqrt(2) * np.sqrt(np.sqrt((kz0**2-kc1**2)**2 + 4 * beta1**2 * k0**4) + kc1**2 - kz0**2)
    kz1 = p1 + 1j * q1
    
    p2 = 1/np.sqrt(2) * np.sqrt(np.sqrt((kz0**2-kc2**2)**2 + 4 * beta2**2 * k0**4) - kc2**2 + kz0**2)
    q2 = 1/np.sqrt(2) * np.sqrt(np.sqrt((kz0**2-kc2**2)**2 + 4 * beta2**2 * k0**4) + kc2**2 - kz0**2)
    kz2 = p2 + 1j * q2
    
    # Reflection coefficients
    r01 = (kz0 - kz1) / (kz0 + kz1)
    r12 = (kz1 - kz2) / (kz1 + kz2)
    
    # Depth grid in film
    z = np.linspace(0, sampthick_angstrom, 1000)
    r01 = r01[:, None]
    r12 = r12[:, None]
    kz1 = kz1[:, None]
    z = z[None, :]
    
    # Calculate electric field intensity
    EE = (1 + r01) * (np.exp(1j * kz1 * z) + r12 * np.exp(1j * kz1 * (2 * sampthick_angstrom - z))) / \
         (1 + r01 * r12 * np.exp(1j * 2 * kz1 * sampthick_angstrom))
    
    if plot:
        xefi_mag = np.abs(EE.T)  # Electric field is complex, get magnitude
        depth = z.flatten()/10
        aois = np.rad2deg(angresamp)

        # Quick plot check 
        cmin = np.quantile(xefi_mag, 0.05)
        cmax = np.quantile(xefi_mag, 1)


        # Plotting
        fig, ax = plt.subplots(figsize=(5,3), dpi=150, tight_layout=True)

        im = ax.imshow(
            xefi_mag, 
            origin='upper', 
            extent=[aois[0],aois[-1],depth[-1],depth[0]], 
            aspect='auto',
            norm=plt.Normalize(cmin,cmax),
            cmap=plt.cm.terrain
        )
        ax.set(xlabel='Incident angle [°]', ylabel= 'Film depth [nm]')
        # ax.yaxis.set_minor_locator(MultipleLocator(10))  # from matplotlib.ticker import MultipleLocator
        # ax.yaxis.set_major_locator(MultipleLocator(20))  # to set specific tick intervals

        # Add colorbar
        cbar = fig.colorbar(im, ax=ax)
        cbar.set_label('XEFI', rotation=270, labelpad=10)
    
    return np.rad2deg(angresamp), z.flatten()/10, EE.T

def calc_critical_angle_table(energies, stoichiometries, densities):
    """
    Calculate and display a table of critical angles for multiple materials and energies.

    Args:
        energies (list): List of X-ray energies in eV.
        stoichiometries (list): List of chemical formulas for different materials.
        densities (list): List of densities corresponding to each material in g/cm³.

    Returns:
        None: Prints a table showing the critical angles for each material and energy.
    """
    all_crit_angles = []
    for i, stoichiometry in enumerate(stoichiometries):
        mat_crit_angles = []
        for energy in energies:
            crit_angle = calc_critical_angle(energy, stoichiometry, densities[i], verbose=False)
            mat_crit_angles.append(crit_angle)
        all_crit_angles.append(mat_crit_angles)
    all_crit_angles = np.asarray(all_crit_angles)
    np.shape(all_crit_angles)
    
    #plot table
    table = PrettyTable()
    table.add_column('Energies (eV)', energies)
    for i, stoichiometry in enumerate(stoichiometries):
        table.add_column(stoichiometry, all_crit_angles[i,:])
        table.float_format = '.4'
    print(f'Critical Angles in degrees')
    print(table)

def calc_yoneda_material(incident_deg, energy, stoichiometry, density, verbose=True):
    """
    Calculate the Yoneda peak position for a given material and incident angle.

    Args:
        incident_deg (float): The incident angle in degrees.
        energy (float): The energy of the X-ray beam in eV.
        stoichiometry (str): The chemical formula of the material.
        density (float): The density of the material in g/cm³.
        verbose (bool, optional): If True, prints the Yoneda peak position. Default is True.

    Returns:
        float: The Yoneda peak position in Å^-1.
    """
    wavelength = 12400/energy # in Å
    crit_deg = calc_critical_angle(energy, stoichiometry, density, verbose=False)
    yoneda_q = 2*np.pi/wavelength*(np.sin(np.deg2rad(incident_deg)) + np.sin(np.deg2rad(crit_deg)))
    if verbose:
        print(f'Yoneda peak for {stoichiometry} at {energy}eV and angle {incident_deg:0.2f}° is {yoneda_q:0.3f}Å^-1')
    return yoneda_q
    
def calc_yoneda_critical(incident_deg, energy, crit_deg, verbose=True):
    """
    Calculate the Yoneda peak position using a specified critical angle.

    Args:
        incident_deg (float): The incident angle in degrees.
        energy (float): The energy of the X-ray beam in eV.
        crit_deg (float): The critical angle in degrees.
        verbose (bool, optional): If True, prints the Yoneda peak position. Default is True.

    Returns:
        float: The Yoneda peak position in Å^-1.
    """
    wavelength = 12400/energy # in Å
    yoneda_q = 2*np.pi/wavelength*(np.sin(np.deg2rad(incident_deg)) + np.sin(np.deg2rad(crit_deg)))
    if verbose:
        print(f'Yoneda peak for critical angle {critical_angle}° and indicent angle {incident_deg:0.2f}° is {yoneda_q:0.3f}Å^-1')
    return yoneda_q

def find_xray_edges(element, verbose=True):
    """
    Find the X-ray absorption edges for a specified element.

    Args:
        element (str): The symbol of the element (e.g., 'Fe' for iron).
        verbose (bool, optional): If True, prints the X-ray edges. Default is True.

    Returns:
        dict: A dictionary containing the edge names and corresponding energies in eV.
    """
    xray_edges = xraydb.xray_edges(element)
    table = PrettyTable()
    table.field_names = ["Edge", "Energy (eV)"]
    # Adding rows to the table
    for edge, data in xray_edges.items():
        table.add_row([edge, data[0]])
    table.float_format = '.0'
    print(f'X-ray Edges for {element}')
    print(table)

    return xray_edges

def calc_mu(energy, stoichiometry, density, verbose=True):
    """
    Calculate the linear attenuation coefficient (mu) for a material at a specific energy.

    Args:
        energy (float): The energy of the X-ray beam in eV.
        stoichiometry (str): The chemical formula of the material.
        density (float): The density of the material in g/cm³.
        verbose (bool, optional): If True, prints the attenuation coefficient. Default is True.

    Returns:
        float: The linear attenuation coefficient (mu) in cm^-1.
    """
    xraydb.add_material('foo', stoichiometry, density)
    mu_val = xraydb.material_mu('foo', energy)
    if verbose:
        print(f'Yoneda peak for {stoichiometry} at {energy:.0f}eV and angle {incident_deg:0.2f}° is {yoneda_q:0.3f}Å^-1')
        
    return mu_val
    
def calc_mu_list(energies, stoichiometry, density):
    """
    Calculate the linear attenuation coefficients (mu) for a material over a range of energies.

    Args:
        energies (list): List of X-ray energies in eV.
        stoichiometry (str): The chemical formula of the material.
        density (float): The density of the material in g/cm³.

    Returns:
        list: A list of linear attenuation coefficients (mu) in cm^-1 for each energy.
    """
    xraydb.add_material('foo', stoichiometry, density)
    mu_vals = xraydb.material_mu('foo', energies)
        
    return mu_vals
    
def calc_mu_grid(energies, stoichiometries, densities, verbose=True, table=True, plot=False):
    """
    Calculate and display the attenuation lengths for multiple materials over a range of energies.

    Args:
        energies (list): List of X-ray energies in eV.
        stoichiometries (list): List of chemical formulas for different materials.
        densities (list): List of densities corresponding to each material in g/cm³.
        verbose (bool, optional): If True, prints the attenuation lengths. Default is True.
        table (bool, optional): If True, displays the attenuation lengths in a table. Default is True.
        plot (bool, optional): If True, plots the attenuation lengths vs. energy. Default is False.

    Returns:
        None: visualizes data in plot or table
    """
    material_mus = []
    for i, stoichiometry in enumerate(stoichiometries):
        mu_vals = calc_mu_list(energies, stoichiometry, densities[i])
        material_mus.append(mu_vals)
    material_mus = np.asarray(material_mus)
    
    #calculate attenuation length (mm)
    material_atts = 10/material_mus
    
    if table:
        #plot table
        table = PrettyTable()
        table.add_column('Energies (eV)', energies)
        for i, stoichiometry in enumerate(stoichiometries):
            table.add_column(stoichiometry, material_atts[i,:])
            table.float_format = '.4'
        print(f'Attenuation lengths in mm')
        print(table)
        
    if plot:
        plt.figure()
        for i, stoichiometry in enumerate(stoichiometries):
            plt.plot(energies, material_atts[i,:], label=stoichiometry)
        plt.xlabel('Energy (eV)')
        plt.ylabel('Attenuation Length (mm)')
        plt.title('Attenuation Length vs Energy')
        plt.legend()
        plt.grid(True)

def calc_transmission(energy, stoichiometry, density, thickness, verbose=True):
    """
    Calculate the transmission of X-rays through a material of specified thickness.

    Args:
        energy (float): The energy of the X-ray beam in eV.
        stoichiometry (str): The chemical formula of the material.
        density (float): The density of the material in g/cm³.
        thickness (float): The thickness of the material in mm.
        verbose (bool, optional): If True, prints the transmission. Default is True.

    Returns:
        float: The transmission fraction (0 to 1) of X-rays through the material.
    """
    mu_val = calc_mu(energy, stoichiometry, density, verbose=False)
    trans_val = np.exp(thicknesses * -0.1*mu_val)
    if verbose:
        print(f'Transmission through {thickness}mm of {stoichiometry} at {energy:.0f}eV is {100*trans_val:0.2f}%')
    return trans_val

def calc_transmission_list(energies, stoichiometry, density, thickness):
    """
    Calculate the transmission of X-rays through a material of specified thickness over a range of energies.

    Args:
        energies (list): List of X-ray energies in eV.
        stoichiometry (str): The chemical formula of the material.
        density (float): The density of the material in g/cm³.
        thickness (float): The thickness of the material in mm.

    Returns:
        list: A list of transmission fractions (0 to 1) for each energy.
    """
    mu_vals = calc_mu_list(energies, stoichiometry, density)
    trans_vals =  np.exp(thickness * -0.1*mu_vals)
        
    return trans_vals

def calc_transmission_grid(energies, stoichiometries, densities, thicknesses, verbose=True, table=True, plot=False):
    """
    Calculate and display the transmission of X-rays through multiple materials and thicknesses over a range of energies.

    Args:
        energies (list): List of X-ray energies in eV.
        stoichiometries (list): List of chemical formulas for different materials.
        densities (list): List of densities corresponding to each material in g/cm³.
        thicknesses (list): List of thicknesses for each material in mm.
        verbose (bool, optional): If True, prints the transmission values. Default is True.
        table (bool, optional): If True, displays the transmission values in a table. Default is True.
        plot (bool, optional): If True, plots the transmission vs. energy. Default is False.

    Returns:
        None: visualizes data in plot or table
    """
    material_trans = []
    for i, stoichiometry in enumerate(stoichiometries):
        trans_vals = calc_transmission_list(energies, stoichiometry, densities[i], thicknesses[i])
        material_trans.append(trans_vals)
    material_trans = np.asarray(material_trans)

    total_trans_vals = []
    for col in range(np.shape(material_trans)[1]):
        total_trans = 1
        for row in range(np.shape(material_trans)[0]):
            total_trans *= material_trans[row,col]
        total_trans_vals.append(total_trans)
    total_trans_vals = np.asarray(total_trans_vals)
    
    col_names = []
    for i, stoichiometry in enumerate(stoichiometries):
        col_name = stoichiometry+' ('+str(thicknesses[i])+' mm)'
        col_names.append(col_name)

    if table:
        #plot table
        table = PrettyTable()
        table.add_column('Energies (eV)', energies)
        for i, stoichiometry in enumerate(stoichiometries):
            table.add_column(col_names[i], 100*material_trans[i,:])
        table.add_column('Total Transmission', 100*total_trans_vals)
        table.float_format = '.2'
        print(f'Transmission values in %')
        print(table)
        
    if plot:
        plt.figure()
        for i, stoichiometry in enumerate(stoichiometries):
            plt.plot(energies, 100*material_trans[i,:], label=f'{stoichiometry}@{thicknesses[i]}mm')
        plt.plot(energies, 100*total_trans_vals, label='Total Transmission', color='black')
        plt.xlabel('Energy (eV)')
        plt.ylabel('Transmission (%)')
        plt.title('X-ray Transmission vs Energy')
        plt.legend()
        plt.grid(True)

def calc_solution_transmission(energies, solvent_stoich, solvent_dens, solute_stoich, solute_dens, mg_per_ml, cap_diam, table=True, plot=False):
    """
    Calculate the transmission of X-rays through a solution in a capillary.

    Args:
        energies (list): List of X-ray energies in eV.
        solvent_stoich (str): The chemical formula of the solvent.
        solvent_dens (float): The density of the solvent in g/cm³.
        solute_stoich (str): The chemical formula of the solute.
        solute_dens (float): The density of the solute in g/cm³.
        mg_per_ml (float): The concentration of the solute in mg/mL.
        cap_diam (float): The diameter of the capillary in mm.
        table (bool, optional): If True, displays the transmission values in a table. Default is True.
        plot (bool, optional): If True, plots the transmission vs. energy. Default is False.

    Returns:
        None: visualizes data in plot or table
    """
    #calculate mu (cm-1)
    solute_vf = (mg_per_ml/1000)/solute_dens
    solute_eff_length = solute_vf*cap_diam
    solv_eff_length = (1-solute_vf)*cap_diam
    solvent_mus = calc_mu_list(energies, solvent_stoich, solvent_dens)
    solute_mus = calc_mu_list(energies, solute_stoich, solute_dens)
    
    #calculate transmission
    solvent_trans = np.exp(solv_eff_length * -0.1*solvent_mus)
    solute_trans = np.exp(solute_eff_length * -0.1*solute_mus)
    total_trans = solvent_trans*solute_trans

    if table:
        #plot table
        table = PrettyTable()
        table.add_column('Energies (eV)', energies)
        table.add_column('Solvent Transmission (%)', solvent_trans*100)
        table.add_column('Solute Transmission (%)', solute_trans*100)
        table.add_column('Total Transmission (%)', total_trans*100)
        table.float_format = '.2'
        print(f'Transmissions through a {cap_diam:.2f}mm capillary with solution of {mg_per_ml:.2f}mg/ml {solute_stoich} dissolved in {solvent_stoich}.')
        print(table)
    if plot:
        plt.figure()
        plt.plot(energies, 100*solvent_trans, label=f'Solvent: {solvent_stoich}')
        plt.plot(energies, 100*solute_trans, label=f'Solute: {solute_stoich}@{mg_per_ml}mg/ml')
        plt.plot(energies, 100*total_trans, label='Total Transmission', color='black')
        plt.xlabel('Energy (eV)')
        plt.ylabel('Transmission (%)')
        plt.title(f'X-ray transmission through {cap_diam}mm capillary')
        plt.legend()
        plt.grid(True)