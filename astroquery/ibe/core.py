# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
IBE
===

API from

 http://irsa.ipac.caltech.edu/ibe/
"""
from __future__ import print_function, division

import os
import warnings

import astropy.units as u
import astropy.coordinates as coord
from astropy.table import Table

from ..exceptions import InvalidQueryError
from ..query import BaseQuery
from ..utils import commons
from . import conf

__all__ = ['Ibe', 'IbeClass']


class IbeClass(BaseQuery):
    URL = conf.server
    MISSION = conf.mission
    DATASET = conf.dataset
    TABLE = conf.table
    TIMEOUT = conf.timeout

    def query_region(
            self, coordinate=None, where=None, mission=None, dataset=None,
            table=None, columns=None, width=None, height=None,
            intersect='OVERLAPS', most_centered=False):
        """
        For certain missions, this function can be used to search for image and
        catalog files based on a point, a box (bounded by great circles) and/or
        an SQL-like ``where`` clause.

        If `coordinates` is specified, then the optional `width` and `height`
        arguments control the width and height of the search box. If neither
        `width` nor `height` are provided, then the search area is a point. If
        only one of `width` or `height` are specified, then the search area is
        a square with that side length centered at the coordinate.

        Parameters
        ----------
        coordinate : str, `astropy.coordinates` object
            Gives the position of the center of the box if performing a box
            search. If it is a string, then it must be a valid argument to
            `astropy.coordinates.SkyCoord`. Required if `where` is absent.
        where : str
            SQL-like query string. Required if `coordinates` is absent.
        mission : str
            The mission to be used (if not the default mission).
        dataset : str
            The dataset to be used (if not the default dataset).
        table : str
            The table to be queried (if not the default table).
        columns : str, list
            A space-separated string or a list of strings of the names of the
            columns to return.
        width : str or `~astropy.units.Quantity` object
            Width of the search box if `coordinates` is present.

            The string must be parsable by `~astropy.coordinates.Angle`. The
            appropriate `~astropy.units.Quantity` object from `astropy.units`
            may also be used.
        height : str, `~astropy.units.Quantity` object
            Height of the search box if `coordinates` is present.

            The string must be parsable by `~astropy.coordinates.Angle`. The
            appropriate `~astropy.units.Quantity` object from `astropy.units`
            may also be used.
        intersect : ``COVERS``, ``ENCLOSED``, ``CENTER``, ``OVERLAPS``
            Spatial relationship between search box and image footprint.

            ``COVERS``: X must completely contain S. Equivalent to ``CENTER``
            and ``OVERLAPS`` if S is a point.

            ``ENCLOSED``: S must completely contain X. If S is a point, the
            query will always return an empty image table.

            ``CENTER``: X must contain the center of S. If S is a point, this
            is equivalent to ``COVERS`` and ``OVERLAPS``.

            ``OVERLAPS``: The intersection of S and X is non-empty. If S is a
            point, this is equivalent to ``CENTER`` and ``COVERS``.
        most_centered : bool
            If True, then only the most centered image is returned.

        Returns
        -------
        table : `~astropy.table.Table`
            A table containing the results of the query
        """
        response = self.query_region_async(
            coordinate=coordinate, where=where, mission=mission,
            dataset=dataset, table=table, columns=columns, width=width,
            height=height, intersect=intersect, most_centered=most_centered)

        # Rause exception, if request failed
        response.raise_for_status()

        return Table.read(response.content, format='ipac')

    def query_region_async(
            self, coordinate=None, where=None, mission=None, dataset=None,
            table=None, columns=None, width=None, height=None,
            intersect='OVERLAPS', most_centered=False):
        """
        For certain missions, this function can be used to search for image and
        catalog files based on a point, a box (bounded by great circles) and/or
        an SQL-like ``where`` clause.

        If `coordinates` is specified, then the optional `width` and `height`
        arguments control the width and height of the search box. If neither
        `width` nor `height` are provided, then the search area is a point. If
        only one of `width` or `height` are specified, then the search area is
        a square with that side length centered at the coordinate.

        Parameters
        ----------
        coordinate : str, `astropy.coordinates` object
            Gives the position of the center of the box if performing a box
            search. If it is a string, then it must be a valid argument to
            `astropy.coordinates.SkyCoord`. Required if `where` is absent.
        where : str
            SQL-like query string. Required if `coordinates` is absent.
        mission : str
            The mission to be used (if not the default mission).
        dataset : str
            The dataset to be used (if not the default dataset).
        table : str
            The table to be queried (if not the default table).
        columns : str, list
            A space-separated string or a list of strings of the names of the
            columns to return.
        width : str or `~astropy.units.Quantity` object
            Width of the search box if `coordinates` is present.

            The string must be parsable by `~astropy.coordinates.Angle`. The
            appropriate `~astropy.units.Quantity` object from `astropy.units`
            may also be used.
        height : str, `~astropy.units.Quantity` object
            Height of the search box if `coordinates` is present.

            The string must be parsable by `~astropy.coordinates.Angle`. The
            appropriate `~astropy.units.Quantity` object from `astropy.units`
            may also be used.
        intersect : ``COVERS``, ``ENCLOSED``, ``CENTER``, ``OVERLAPS``
            Spatial relationship between search box and image footprint.

            ``COVERS``: X must completely contain S. Equivalent to ``CENTER``
            and ``OVERLAPS`` if S is a point.

            ``ENCLOSED``: S must completely contain X. If S is a point, the
            query will always return an empty image table.

            ``CENTER``: X must contain the center of S. If S is a point, this
            is equivalent to ``COVERS`` and ``OVERLAPS``.

            ``OVERLAPS``: The intersection of S and X is non-empty. If S is a
            point, this is equivalent to ``CENTER`` and ``COVERS``.
        most_centered : bool
            If True, then only the most centered image is returned.

        Returns
        -------
        response : `requests.Response`
            The HTTP response returned from the service
        """

        if coordinate is None and where is None:
            raise InvalidQueryError(
                'At least one of `coordinate` or `where` is required')

        intersect = intersect.upper()
        if intersect not in ('COVERS', 'ENCLOSED', 'CENTER', 'OVERLAPS'):
            raise InvalidQueryError(
                "Invalid value for `intersects` " +
                "(must be 'COVERS', 'ENCLOSED', 'CENTER', or 'OVERLAPS')" )

        args = {
            'INTERSECT': intersect
        }

        # Note: in IBE, if 'mcen' argument is present, it is true.
        # If absent, it is false.
        if most_centered:
            args['mcen'] = '1'

        if coordinate is not None:
            c = commons.parse_coordinates(coordinate).transform_to(coord.ICRS)
            args['POS'] = '{0},{1}'.format(c.ra.deg, c.dec.deg)
            if width and height:
                args['SIZE'] = '{0},{1}'.format(
                    commons.parse_radius(width).value,
                    commons.parse_radius(height).value)
            elif width or height:
                args['SIZE'] = str(commons.parse_radius(width or height).value)

        if where:
            args['where'] = where

        if columns:
            if isinstance(columns, basestring):
                columns = columns.split()
            args['columns'] = ','.join(columns)

        url = os.path.join(
            self.URL, 'search',
            mission or self.MISSION,
            dataset or self.DATASET,
            table or self.TABLE)

        return commons.send_request(url, args, self.TIMEOUT, request_type='GET')

Ibe = IbeClass()
