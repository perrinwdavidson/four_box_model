class FourBox:
    """
    purpose :: to solve a coupled ocean-atmosphere model, as presented in Sarmiento and Toggweiler (1984, 1985)
    author :: perrin w. davidson
    contact :: perrinwdavidson@gmail.com
    date :: 5.12.2022
    """
    # import libraries ::
    import numpy as np

    # define initial values, taken from Emerson and Hamme (2022) ::
    V0 = 1.292E18  # ocean volume [m+3]
    S0 = 3.49E14  # ocean surface area [m+2]
    Sl_perc = 0.85  # percentage of surface ocean low latitude [%]
    Sl = S0 * Sl_perc  # surface area of low latitude box [m+2]
    Sh_perc = 0.15  # percentage of surface ocean high latitude [%]
    Sh = S0 * Sh_perc  # surface area of high latitude box [m+2]
    Depth_l = 100  # depth of low latitude box [m]
    Vl = Sl * Depth_l  # volume of low latitude box [m+3]
    Depth_h = 250  # depth of high latitude box [m]
    Vh = Sh * Depth_h  # volume of high latitude box [m+3]
    Vd = V0 - Vh - Vl  # volume of deep ocean [m+3]
    Ma = 1.773E20  # molar volume of the atmosphere [mol atm-1]
    Temp_l = 21.5  # temperature of low latitude box [C]
    Temp_h = 2.5  # temperature of high latitude box [C]
    Sal = 34.7  # salinity of surface ocean []
    K_H_l = 31.9  # henry's law coefficient for low latitude co2 [mol m-3 atm-1]
    K_H_h = 58.8  # henry's law coefficient for high latitude co2 [mol m-3 atm-1]
    Pd = 2.15E-3  # phosphate concentration of the deep box [mol m-3]
    Sum_C = 3.02E18  # total C in the ocean (dic) and atmosphere (fco2) [mol]
    Sum_A = 3.14E18  # total alkalinity in the ocean [eq]
    T = 20E6  # thermohaline circulation [m+3 s-1]
    r_ap = 44  # ratio of alkalinity to phosphorous in biological flux term
    r_cp = 136  # ratio of carbon to phosphorous in biological flux term
    k = 3  # gas exchange mas transfer coefficient [m d-1]
    beta_l = -784.16  # linear relationship of co2, alk, and dic: (Ai - Ci) = (beta_i * fco2_i) + gamma_i
    gamma_l = 0.56759  # linear relationship of co2, alk, and dic: (Ai - Ci) = (beta_i * fco2_i) + gamma_i
    beta_h = -620.04  # linear relationship of co2, alk, and dic: (Ai - Ci) = (beta_i * fco2_i) + gamma_i
    gamma_h = 0.37159  # linear relationship of co2, alk, and dic: (Ai - Ci) = (beta_i * fco2_i) + gamma_i
    scale_fact = 1E15  # scalar factor to remove problems where some rows of carbon matrix contain large numbers

    # set values to be solved ::
    Pl = np.nan
    Ph = np.nan
    Al = np.nan
    Ah = np.nan
    Ad = np.nan
    Cl = np.nan
    Ch = np.nan
    Cd = np.nan
    fCO2_l = np.nan
    fCO2_h = np.nan
    fCO2_a = np.nan

    # solve phosphate steady state mass balance equations ::
    def solve_po4(self, j, f_hd):
        """
        purpose :: to solve steady state mass balance equations for phosphate
        inputs ::
            [1] j - biological carbon flux pump flux of phosphate [mol P s-1]
            [2] f_hd - high latitude-deep box exchange rate [m+3 s-1]
        outputs ::
            [1] x_po4 - steady state concentrations of phosphate in the form [high lat, low lat]
        """
        # import libraries ::
        import numpy as np

        # define matrices ::
        a_po4 = np.array([[(-f_hd - self.T), self.T],
                         [0, -self.T]])
        b_po4 = np.array([[(j[0] * self.Sh) - (f_hd * self.Pd)],
                          [(j[1] * self.Sl) - (self.T * self.Pd)]])

        # solve for concentrations ::
        x_po4 = np.linalg.solve(a_po4, b_po4)

        # set values ::
        self.Ph, self.Pl = x_po4

        # return concentrations ::
        return x_po4

    # solve alkalinity in steady-steady state ::
    def solve_alk(self, j, f_hd):
        """
        purpose :: to solve steady state mass balance equations for alkalinity
        inputs ::
            [1] j - biological carbon flux pump flux of phosphate [mol P s-1]
            [2] f_hd - high latitude-deep box exchange rate [m+3 s-1]
        outputs ::
            [1] x_alk - steady state concentrations of alkalinity in the form [low lat, high lat, deep]
        """
        # import libraries ::
        import numpy as np

        # define matrices ::
        a_alk = np.array([[self.T, (-self.T - f_hd), f_hd],
                          [-self.T, 0, self.T],
                          [self.Vl, self.Vh, self.Vd]])
        b_alk = np.array([[j[0] * self.Sh * self.r_ap],
                          [j[1] * self.Sl * self.r_ap],
                          [self.Sum_A]])

        # solve for concentrations ::
        x_alk = np.linalg.solve(a_alk, b_alk)

        # set values ::
        self.Al, self.Ah, self.Ad = x_alk

        # return concentrations ::
        return x_alk

    # solve co2 in steady-steady state ::
    def solve_co2(self, j, f_hd):
        """
        purpose :: to solve steady state mass balance equations for co2 in ocean and atmosphere with linear relation
        inputs ::
            [1] j - biological carbon flux pump flux of phosphate [mol P s-1]
            [2] f_hd - high latitude-deep box exchange rate [m+3 s-1]
        outputs ::
            [1] x_co2 - steady state concentrations of co2 in the form [low lat, high lat, deep] * [ocean, atmosphere]
        """
        # import libraries ::
        import numpy as np

        # define matrices ::
        a_co2 = np.array([[self.T, (-self.T - f_hd), f_hd, 0, (-self.k * self.K_H_h * self.Sh), (self.k * self.K_H_h * self.Sh)],
                          [-self.T, 0, self.T, (-self.k * self.K_H_l * self.Sl), 0, (self.k * self.K_H_l * self.Sl)],
                          [0, 0, 0, (self.k * self.K_H_l * self.Sl), (self.k * self.K_H_h * self.Sh), (-self.k * ((self.K_H_h * self.Sh) + (self.K_H_l * self.Sl)))],
                          [self.Vl, self.Vh, self.Vd, 0, 0, self.Ma],
                          [self.scale_fact, 0, 0, (self.beta_l * self.scale_fact), 0, 0],
                          [0, self.scale_fact, 0, 0, (self.beta_h * self.scale_fact), 0]])
        b_co2 = np.array([[(j[0] * self.Sh * self.r_cp)],
                          [(j[1] * self.Sl * self.r_cp)],
                          [0],
                          [self.Sum_C],
                          [(self.Al[0] - self.gamma_l) * self.scale_fact],
                          [(self.Ah[0] - self.gamma_h) * self.scale_fact]])

        # solve for concentrations ::
        x_co2 = np.linalg.solve(a_co2, b_co2)

        # set values ::
        self.Cl, self.Ch, self.Cd, self.fCO2_l, self.fCO2_h, self.fCO2_a = x_co2

        # return concentrations ::
        return x_co2

    # run experiments ::
    def contour_data(self, j_range, f_hd_range):
        """
        purpose :: to solve steady state mass balance equations for alkalinity
        inputs ::
            [1] j_range - range of biological carbon flux pump flux of phosphate [mol P s-1]
            [2] f_hd_range - range of high latitude-deep box exchange rate [m+3 s-1]
        outputs ::
            [0] plot - plot of contours of fco2, dip, and efficiency over ranges specified
            [1] data_countours - data from contouring in the form [fco2a, dip_h, efficiency]
        """
        # import libraries ::
        import numpy as np
        import matplotlib.pyplot as plt

        # set conversions ::
        y2s = 60 ** 2 * 24 * 365
        c2p = 106

        # set j low ::
        jl = 0.15 / y2s / c2p

        # preallocate array ::
        data_contours = np.empty((len(j_range), len(f_hd_range), 3))
        data_contours[:] = np.nan

        # run experiment ::
        for ij in range(len(j_range)):
            for if_hd in range(len(f_hd_range)):

                # calculate po4_h ::
                ix_po4 = self.solve_po4([j_range[ij], jl], f_hd_range[if_hd])
                val_po4 = ix_po4[0] * 1E6 / 1028
                if val_po4 <= 0:
                    val_po4 = np.nan
                data_contours[ij, if_hd, 1] = val_po4

                # calculate co2_a ::
                ix_co2 = self.solve_co2([j_range[ij], jl], f_hd_range[if_hd])
                val_co2 = ix_co2[5] * 1E6
                if np.isnan(val_po4):
                    val_co2 = np.nan
                data_contours[ij, if_hd, 0] = val_co2

        # calculate efficiencies ::
        data_contours[:, :, 2] = 1 - (data_contours[:, :, 1] / (self.Pd * 1E6 / 1028))

        # make plotting arrays ::
        xplot, yplot = np.meshgrid(f_hd_range / 1E6, j_range * y2s * c2p)

        # set plotting settings ::
        plt.rc('text', usetex=True)
        plt.rcParams.update({'font.size': 16})
        plt.rc('font', family='serif')
        plt.rc('text.latex', preamble=r'\usepackage{wasysym}')

        # plot ::
        fig, ax = plt.subplots(3, 1, figsize=(10, 24))
        fig1 = ax[0].contourf(xplot, yplot, data_contours[:, :, 0])
        ax[0].set_xlabel("$f_{hd}$ (Sverdrups)")
        ax[0].set_ylabel("$J_{h}$ (mol C m$^{-2}$ y$^{-1})$")
        plt.colorbar(fig1, ax=ax[0], label="$f$CO$_2^a$ ($\mu$atm)")
        fig2 = ax[1].contourf(xplot, yplot, data_contours[:, :, 1])
        ax[1].set_xlabel("$f_{hd}$ (Sverdrups)")
        ax[1].set_ylabel("$J_{h}$ (mol C m$^{-2}$ y$^{-1})$")
        plt.colorbar(fig2, ax=ax[1], label="DIP ($\mu$mol kg$^{-1}$)")
        fig3 = ax[2].contourf(xplot, yplot, data_contours[:, :, 2], vmin=0, vmax=1)
        ax[2].set_xlabel("$f_{hd}$ (Sverdrups)")
        ax[2].set_ylabel("$J_{h}$ (mol C m$^{-2}$ y$^{-1})$")
        plt.colorbar(fig3, ax=ax[2], label="Efficiency ($\%$)", ticks=[0, 0.2, 0.4, 0.6, 0.8, 1.0])
        plt.savefig("plots/experiment1.png", dpi=300)
        plt.show()

        # return data ::
        return data_contours

# end class
