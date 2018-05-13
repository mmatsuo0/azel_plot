#!/usr/bin/env python
# Python 2.7.14

import argparse
import numpy
import pandas
import astropy.coordinates
import astropy.units as u
import matplotlib.pyplot
import os

site = 'Nobeyama'
day_default = '2018-06-06'

matplotlib.rcParams['lines.markersize'] = 1
matplotlib.rcParams['font.family'] = 'Times New Roman'
matplotlib.rcParams['axes.grid'] = True
matplotlib.rcParams['grid.linestyle'] = ':'
matplotlib.rcParams['scatter.marker'] = '.'


class LST:
    def __init__(self, file_name, day):
        self.objects_file = file_name
        if day:
            self.day = day
            self.sun = True
        else:
            self.day = day_default
            self.sun = False

    def get_params(self, s):
        if s == 'Nobeyama':
            lat = 35.9446944444 * u.deg
            lon = 138.4725555556 * u.deg
            height = 1350. * u.m

        location = astropy.coordinates.EarthLocation(lat=lat, lon=lon,
                                                     height=height)
        utc_t = astropy.time.Time(self.day + ' 00:00:00')
        lst_t = utc_t.sidereal_time('mean', longitude=lon)
        dlst = (24. - 3. / 60. - lst_t.hour) * u.hour
        utc1 = utc_t + dlst
        oneday = numpy.linspace(5. / 60, 24 - 5. / 60., 1000) * u.hour
        self.utc_1day = utc1 + oneday
        lst_1day = self.utc_1day.sidereal_time('mean', longitude=lon)
        self.lst_1day = lst_1day.hour

        self.altazframe = astropy.coordinates.AltAz(obstime=self.utc_1day, location=location)

    def plot(self):
        objects = pandas.read_csv(self.objects_file, comment='#')
        file_base = os.path.splitext(self.objects_file)[0]

        fig = matplotlib.pyplot.figure(figsize=(10, 8))
        ax1 = fig.add_subplot(211)
        ax2 = fig.add_subplot(212)

        if self.sun:
            sunaltaz = astropy.coordinates.get_sun(self.utc_1day).transform_to(self.altazframe)
            LST = self.lst_1day[numpy.where(sunaltaz.alt.deg > 0)]
            AZ = sunaltaz.az.deg[numpy.where(sunaltaz.alt.deg > 0)]
            ax1.scatter(LST, AZ, c='r')
            ax2.plot(self.lst_1day, sunaltaz.alt.deg, c='r', label='Sun')
            fig.text(0.2, 0.9, self.day)

        for i in range(len(objects)):
            c = astropy.coordinates.SkyCoord(objects['X'].iloc[i],
                                             objects['Y'].iloc[i],
                                             frame=objects['frame'].iloc[i])
            altaz = c.transform_to(self.altazframe)
            LST = self.lst_1day[numpy.where(altaz.alt.deg > 0)]
            AZ = altaz.az.deg[numpy.where(altaz.alt.deg > 0)]
            ax1.scatter(LST, AZ, c=objects['color'].iloc[i])
            ax2.plot(self.lst_1day, altaz.alt.deg, c=objects['color'].iloc[i],
                     label=objects['OBJECT'].iloc[i])

        ax1.xaxis.set_major_locator(matplotlib.ticker.MultipleLocator(1))
        ax1.yaxis.set_major_locator(matplotlib.ticker.MultipleLocator(30))
        ax1.set_xlim(0, 24)
        ax1.set_ylim(0, 360)
        ax1.set_ylabel('AZ (degree)')
        ax2.xaxis.set_major_locator(matplotlib.ticker.MultipleLocator(1))
        ax2.set_xlim(0, 24)
        ax2.set_ylim(0, 90)
        ax2.axhline(12, c='r', ls='--')
        ax2.axhline(80, c='r', ls='--')
        ax2.set_xlabel('LST (h)')
        ax2.set_ylabel('EL (degree)')

        ax2.legend(loc='lower left', ncol=5, bbox_to_anchor=(0.3, 2.2), fontsize=8)

        if self.sun:
            fig.savefig('{}_LST_{}.png'.format(file_base, self.day))
        else:
            fig.savefig('{}_LST.png'.format(file_base))


class JST:
    def __init__(self, file_name, day):
        self.objects_file = file_name
        if day:
            self.day = day
        else:
            self.day = day_default

    def get_params(self, s):
        if s == 'Nobeyama':
            lat = 35.9446944444 * u.deg
            lon = 138.4725555556 * u.deg
            height = 1350. * u.m
        location = astropy.coordinates.EarthLocation(lat=lat, lon=lon, height=height)

        utcoffset = 9 * u.hour
        jst1 = astropy.time.Time(self.day + ' 00:00:00')
        utc1 = jst1 - utcoffset
        oneday = numpy.linspace(0, 24, 1000) * u.hour
        self.utc_1day = utc1 + oneday
        self.jst_1day = jst1 + oneday

        self.altazframe = astropy.coordinates.AltAz(obstime=self.utc_1day, location=location)

    def plot(self):
        objects = pandas.read_csv(self.objects_file, comment='#')
        file_base = os.path.splitext(self.objects_file)[0]

        fig = matplotlib.pyplot.figure(figsize=(10, 8))
        ax1 = fig.add_subplot(211)
        ax2 = fig.add_subplot(212)

        sunaltaz = astropy.coordinates.get_sun(self.utc_1day).transform_to(self.altazframe)
        JST = self.jst_1day[numpy.where(sunaltaz.alt.deg > 0)]
        AZ = sunaltaz.az.deg[numpy.where(sunaltaz.alt.deg > 0)]
        ax1.scatter(JST.datetime, AZ, c='r')
        ax2.plot(self.jst_1day.datetime, sunaltaz.alt.deg, c='r', label='Sun')

        for i in range(len(objects)):
            c = astropy.coordinates.SkyCoord(objects['X'].iloc[i],
                                             objects['Y'].iloc[i],
                                             frame=objects['frame'].iloc[i])
            altaz = c.transform_to(self.altazframe)
            JST = self.jst_1day[numpy.where(altaz.alt.deg > 0)]
            AZ = altaz.az.deg[numpy.where(altaz.alt.deg > 0)]
            ax1.scatter(JST.datetime, AZ, c=objects['color'].iloc[i])
            ax2.plot(self.jst_1day.datetime, altaz.alt.deg, c=objects['color'].iloc[i],
                     label=objects['OBJECT'].iloc[i])

        ax1.xaxis.set_major_locator(matplotlib.dates.HourLocator())
        ax1.xaxis.set_major_formatter(matplotlib.dates.DateFormatter('%H'))
        ax1.yaxis.set_major_locator(matplotlib.ticker.MultipleLocator(30))
        ax1.set_xlim(self.jst_1day.datetime[0], self.jst_1day.datetime[-1])
        ax1.set_ylim(0, 360)
        ax1.set_ylabel('AZ (degree)')
        ax2.xaxis.set_major_locator(matplotlib.dates.HourLocator())
        ax2.xaxis.set_major_formatter(matplotlib.dates.DateFormatter('%H'))
        ax2.set_xlim(self.jst_1day.datetime[0], self.jst_1day.datetime[-1])
        ax2.set_ylim(0, 90)
        ax2.axhline(12, c='r', ls='--')
        ax2.axhline(80, c='r', ls='--')
        ax2.set_xlabel('JST (h)')
        ax2.set_ylabel('EL (degree)')

        ax2.legend(loc='lower left', ncol=5, bbox_to_anchor=(0.3, 2.2), fontsize=8)
        fig.text(0.2, 0.9, self.day)

        fig.savefig('{}_JST_{}.png'.format(file_base, self.day))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('file', nargs='?', default='objects.txt')
    parser.add_argument('-d', nargs='?')
    parser.add_argument('--jst', nargs='?', default='LST', const='JST')
    args = parser.parse_args()

    if args.jst == 'LST':
        print 'LST'
        lst = LST(args.file, args.d)
        lst.get_params(site)
        lst.plot()
    elif args.jst == 'JST':
        print 'JST'
        jst = JST(args.file, args.d)
        jst.get_params(site)
        jst.plot()
