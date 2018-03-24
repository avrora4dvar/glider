# glider


## Figure file combinations by figure number

1. Matlab\surface\cal_layer.m - Python\plt1_map.py
2. Matlab\wind\plt_wind_point.m Python\plt1_wind.m
3. Python\plt1_map.py
4. Matlab\article_glider\cal_cov_z.m - Python\cov_z.py
5. Matlab\surface\cal_layer.m  - Python\plt1_vort_plot.py
6. Matlab\surface\cal_layer.m - Python\plt1_temp_plot.py
7. Matlab\cal_cross_alongglider.m - Python\plt1_along_glider.py
8. Matlab\cross\cal_crosslon.m - Python\plt1_cross_rho_plot.py
9. Matlab\article_glider\cal_surface_rho.m - Python\plt1_variance_rho.py
10. Matlab\article_glider\cal_surface_rho.m - Python\plt1_rho_surface_free_plot.py
11. Matlab\article_glider\plt_TSz.m - Python\plt1_zprofile_plot.py
12. Matlab\article_glider\plt_TSz.m - Python\plt1_zprofile_plot.py
13. Matlab\article_glider\cal_stat_long.m - Python\plt1_long_sample.py
14. Matlab\article_glider\cal_rms_glider_variations.m - Python\plt1_glider_variation.py
15. Matlab\cov\cal_ritz_vectors.m - Matlab\article_glider\cal_ritz_vectors.m - Python\plt1_array_modes.py
16. Matlab\article_glider\cal_pot_energy.m - Python\plt1_Epot_plot.py
17. Matlab\article_glider\cal_pot_energy.m - Python\plt1_Epot_plot.py

## Caveats

- For lon,lat-figures that sample the hourly output line 22 in Matlab/get_roms_layer needs to be set to 'true'. 
For figures that use daily-averaged output it must be 'false'
- For lon,z-figures that sample the hourly output line 64 in Matlab/get_roms_cross_zonal needs to be set to 'ocean_his'. 
For figures that use daily-averaged output it must be 'ocean_avg'
