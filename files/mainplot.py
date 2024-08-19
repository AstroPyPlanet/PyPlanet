import numpy as np
import matplotlib.pyplot as plt
import io
import dearpygui.dearpygui as dpg
from . import global_vars
from astropy.io import fits
from astropy.wcs import WCS
from PIL import Image
from astropy.nddata import Cutout2D
from astropy.visualization import AsinhStretch, ImageNormalize
from scipy.interpolate import interp2d

button_visible = False

def show_fits_image(selected_file):
    '''
    Displays a FITS image in a window using Dear PyGui.

    Parameters:
    selected_file (str): Path to the FITS file to be displayed.

    The function reads the FITS file, extracts image data, and converts it into a texture to be shown in a Dear PyGui window.
    Additionally, it allows the user to select coordinates on the image and save those coordinates.
    '''
    hdu = fits.open(selected_file)
    data = 1e3 * np.squeeze(hdu[0].data)
    header = hdu[0].header

    w = WCS(header)
    fig, ax = plt.subplots(figsize=(10, 10))
    im = ax.imshow(data, cmap='plasma', origin='upper', aspect='equal', interpolation='nearest')
    ax.axis('off')

    plt.subplots_adjust(left=0, right=1, top=1, bottom=0)
    plt.margins(0, 0)
    plt.gca().xaxis.set_major_locator(plt.NullLocator())
    plt.gca().yaxis.set_major_locator(plt.NullLocator())

    buf = io.BytesIO()
    fig.savefig(buf, format='png')
    buf.seek(0)
    plt.close(fig)

    image = Image.open(buf)
    width, height = image.size
    image_data = np.array(image).astype(np.float32) / 255.0
    image_data = image_data.flatten()

    with dpg.texture_registry(show=False):
        dpg.add_static_texture(width, height, image_data, tag="plot_texture")

    with dpg.window(label="Disk Window", no_move=True, no_title_bar=True, no_resize=True) as disk_window:
        dpg.add_spacer(height=5)
        dpg.add_text("   Click the center of the disk with the middle mouse button to select the coordinates")
        dpg.add_spacer(height=3)
        dpg.add_text("     Once you have entered the data, click the button to generate the deprojection")

        dpg.add_spacer(height=5)
        with dpg.group(horizontal=True):
            dpg.add_spacer(width=180)
            dpg.add_spacer(height=20)
            dpg.add_button(label="   Generate deprojection   ", callback=lambda: deproject_image(selected_file=selected_file), tag="generate_button", show=False)
        
        dpg.add_spacer(height=5)
        dpg.add_text("    Status: Waiting for action...", tag="status_text", color=(255, 255, 0))
        dpg.add_spacer(height=5)

    dpg.set_item_pos(disk_window, [8, 320])
    dpg.set_item_width(disk_window, 650)
    dpg.set_item_height(disk_window, 670)

    with dpg.plot(label="Disk without deprojection", width=500, height=500, pos=(70,150), parent=disk_window, tag="main_plot"):
        dpg.add_plot_legend()
        x_axis = dpg.add_plot_axis(dpg.mvXAxis, label="X-Axis")
        y_axis = dpg.add_plot_axis(dpg.mvYAxis, label="Y-Axis")
        dpg.add_image_series("plot_texture", [0, height], [width, 0], parent=y_axis)

    with dpg.handler_registry():
        dpg.add_mouse_click_handler(callback=lambda sender, app_data: on_plot_click(sender, app_data))

def on_plot_click(sender, app_data):
    '''
    Handles mouse click events on the plot.

    Parameters:
    sender (int): The ID of the item that triggered the event.
    app_data (int): Event data, including which mouse button was pressed.

    If the middle mouse button is clicked, it shows a popup window to confirm the selected coordinates.
    '''
    plot_x, plot_y = dpg.get_plot_mouse_pos()

    if app_data == dpg.mvMouseButton_Middle:
        if dpg.does_item_exist("popup_window"):
            dpg.delete_item("popup_window")

        with dpg.window(label="Choose the center", modal=True, pos=(190, 650), tag="popup_window", no_resize=True):
            dpg.add_text(f"Selected coordinates: X: {plot_x:.2f}, Y: {plot_y:.2f}")

            with dpg.group(horizontal=True):
                dpg.add_button(label="Accept", callback=lambda: save_coordinates(plot_x, plot_y))
                dpg.add_button(label="Cancel", callback=lambda: dpg.delete_item("popup_window"))
            dpg.add_spacer(height=2)
            dpg.add_text("", id="saved_coords", color=(160, 230, 230))

def save_coordinates(x, y):
    '''
    Saves the selected coordinates to a global variable.

    Parameters:
    x (float): The selected X coordinate.
    y (float): The selected Y coordinate.

    Updates the global variable with the selected coordinates and displays a confirmation message.
    '''
    global button_visible
    global_vars.center_coords_value = [x, y]
    dpg.set_value("saved_coords", "Coordinates saved")

    
    button_visible = True
    dpg.show_item("generate_button")

def deproject_image(selected_file):
    '''
    Performs a deprojection of the FITS image and displays the result in a new window.

    Parameters:
    selected_file (str): Path to the FITS file to be processed.

    The function reads the FITS file, performs a deprojection based on the selected coordinates and other global variables,
    and finally displays the deprojected image in a new window.
    '''
    dpg.set_value("status_text", "    Status: Loading...")
    center_x_1000 = global_vars.center_coords_value[0]
    center_y_1000 = global_vars.center_coords_value[1]
    
    hdu = fits.open(selected_file)
    data = 1e3 * np.squeeze(hdu[0].data)
    header = hdu[0].header
    
    original_size_x = header['NAXIS1']
    original_size_y = header['NAXIS2']

    adjusted_size = 1000
    center_x = center_x_1000 * (original_size_x / adjusted_size)
    center_y = center_y_1000 * (original_size_y / adjusted_size)

    w = WCS(header)
    w_ra_dec = w.sub(2)

    cutout_size = (1200, 1200)

    cutout_center_x = center_x
    cutout_center_y = center_y

    cutout = Cutout2D(data, (cutout_center_x, cutout_center_y), cutout_size, wcs=w_ra_dec)
    image_cutout = cutout.data

    subhd = cutout.wcs.to_header()
    ny, nx = image_cutout.shape
    RAo = 3600 * subhd['CDELT1'] * (np.arange(nx) - (subhd['CRPIX1'] - 1))
    DECo = 3600 * subhd['CDELT2'] * (np.arange(ny) - (subhd['CRPIX2'] - 1))

    center_x_header = header['CRPIX1']
    center_y_header = header['CRPIX2']

    distance_x = center_x - center_x_header
    distance_y = center_y - center_y_header

    scale_factor = 1.5 / 500
    real_x = scale_factor * distance_x * -1
    real_y = scale_factor * distance_y

    offRA, offDEC = real_x, real_y
    RAo_shift, DECo_shift = RAo - offRA, DECo - offDEC

    beam_maj, beam_min, beam_PA = 3600 * header['BMAJ'], 3600 * header['BMIN'], header['BPA']
    beam_area = (np.pi * beam_maj * beam_min / (4 * np.log(2))) / (3600 * 180 / np.pi)**2

    vmin, vmax = 0.0, 4.0
    norm = ImageNormalize(vmin=vmin, vmax=vmax, stretch=AsinhStretch(0.03))

    incl = global_vars.incl_value
    PA = global_vars.pa_value
    inclr, PAr = np.radians(incl), np.radians(PA)
    rout = global_vars.rout

    f = interp2d(RAo_shift, DECo_shift, image_cutout, fill_value=0)
    global_vars.f = f

    disk_delt = 3600 * subhd['CDELT2']
    global_vars.disk_delt = disk_delt
    nxd = nyd = 1001
    i0 = nxd // 2
    j0 = nyd // 2
    xdisk = (np.arange(nxd) - i0) * (-disk_delt)
    global_vars.xdisk = xdisk
    ydisk = (np.arange(nyd) - j0) * disk_delt
    global_vars.ydisk = ydisk
    image_deproj = np.zeros((nyd, nxd))

    for j in range(nyd):
        for i in range(nxd):
            dra = np.cos(PAr) * xdisk[i] * np.cos(inclr) + np.sin(PAr) * ydisk[j]
            ddec = -np.sin(PAr) * xdisk[i] * np.cos(inclr) + np.cos(PAr) * ydisk[j]
            image_deproj[j, i] = f(dra, ddec)
        
    global_vars.image_deproj = image_deproj

    deproj_bounds = (xdisk.min(), xdisk.max(), ydisk.min(), ydisk.max())

    fig, ax = plt.subplots(figsize=(10, 8))
    im = ax.imshow(image_deproj, origin='upper', cmap='magma', extent=deproj_bounds, norm=norm, aspect='equal')
    ax.axis('off')

    plt.subplots_adjust(left=0, right=1, top=1, bottom=0)
    plt.margins(0, 0)

    buf = io.BytesIO()
    fig.savefig(buf, format='png', bbox_inches='tight', pad_inches=0)
    buf.seek(0)
    plt.close(fig)

    image = Image.open(buf)
    width, height = image.size
    image_data = np.array(image).astype(np.float32) / 255.0

    with dpg.texture_registry(show=False):
        dpg.add_static_texture(width, height, image_data, tag="deproj_texture")

    with dpg.window(label="Deprojected Image", no_move=True, no_title_bar=True, no_resize=True) as deproj_window:
        
        with dpg.plot(label="Disk with the Deprojection", width=-1, height=-1, parent=deproj_window, tag="deproj_plot"):
            dpg.add_plot_legend()
            x_axis = dpg.add_plot_axis(dpg.mvXAxis, label="X-Axis")
            y_axis = dpg.add_plot_axis(dpg.mvYAxis, label="Y-Axis")
            dpg.add_image_series("deproj_texture", [0, height], [width, 0], parent=y_axis)

    dpg.set_item_pos(deproj_window, [684, 393])
    dpg.set_item_width(deproj_window, 600)
    dpg.set_item_height(deproj_window, 600)

    global_vars.header = header
    global_vars.subhd = subhd
    global_vars.inclr = inclr
    global_vars.PAr = PAr

    f = interp2d(global_vars.xdisk, global_vars.ydisk, global_vars.image_deproj, kind='cubic', fill_value=0)

    dr = global_vars.disk_delt
    r = np.arange(dr, 2*global_vars.rout, dr)
    theta = np.linspace(-180, 180, 181)
    image_polar = np.zeros((theta.size, r.size))

    for j in range(theta.size):
        for i in range(r.size):
            xd = r[i] * np.cos(np.radians(theta[j]))
            yd = r[i] * np.sin(np.radians(theta[j]))
            image_polar[j,i] = f(xd, yd)

    beam_maj, beam_min, beam_PA = 3600 * global_vars.header['BMAJ'], 3600 * global_vars.header['BMIN'], global_vars.header['BPA']
    beam_area = (np.pi * beam_maj * beam_min / (4 * np.log(2))) / (3600 * 180 / np.pi)**2

    freq = global_vars.header['CRVAL3']

    kB_, c_ = 1.38064852e-16, 2.99792e10
    Tb = c_**2 * 1e-26 * np.mean(image_polar, axis=0) / beam_area / (2 * kB_ * freq**2)
    Tb_err = c_**2 * 1e-26 * np.std(image_polar, axis=0) / np.sqrt(image_polar.shape[1]) / beam_area / (2 * kB_ * freq**2)

    with dpg.window(label="Radial Profile", tag="window_tb", no_move=True, no_resize=True, no_title_bar=True, pos=(684, 25)):

        with dpg.plot(label="Radial Profile", height=340, width=584, tag="plot_tb"):
            x_axis = dpg.add_plot_axis(dpg.mvXAxis, label="r [arcsec]", tag="x_axis_tb")
            dpg.set_axis_limits(x_axis, 0, global_vars.rout)
            y_axis = dpg.add_plot_axis(dpg.mvYAxis, label="T_b [K]", tag="y_axis_tb", log_scale=True)
            dpg.set_axis_limits(y_axis, 0.3, 100)
            dpg.add_line_series(r, Tb, label="Brightness Temperature", parent=y_axis, tag="line_series_tb")
        
    dpg.set_value("status_text", "    Deprojection complete")
