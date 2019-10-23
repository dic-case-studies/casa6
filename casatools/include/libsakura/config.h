/*
 * 
 * Copyright (C) 2013-2016
 * National Astronomical Observatory of Japan
 * 2-21-1, Osawa, Mitaka, Tokyo, 181-8588, Japan.
 * 
 * This file is part of Sakura.
 * 
 * Sakura is free software: you can redistribute it and/or modify it under 
 * the terms of the GNU Lesser General Public License as published by the 
 * Free Software Foundation, either version 3 of the License, or (at your 
 * option) any later version.
 * 
 * Sakura is distributed in the hope that it will be useful, but WITHOUT 
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or 
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public 
 * License for more details.
 * 
 * You should have received a copy of the GNU Lesser General Public License 
 * along with Sakura.  If not, see <http://www.gnu.org/licenses/>.
 * 
 */
#ifndef LIBSAKURA_LIBSAKURA_CONFIG_H_
#define LIBSAKURA_LIBSAKURA_CONFIG_H_

#define LIBSAKURA_VERSION_MAJOR 3
#define LIBSAKURA_VERSION_MINOR 0
#define LIBSAKURA_VERSION_STRING "3.0"
#ifndef LIBSAKURA_PREFIX
# define LIBSAKURA_PREFIX	sakura
# define LIBSAKURA_PREFIX_STRING	"sakura"
#endif
#define LIBSAKURA_CONCAT(x, y)	x ## _ ## y
#define LIBSAKURA_SYMBOL(x)	LIBSAKURA_CONCAT(sakura, x)

#define LIBSAKURA_HAS_LOG4CXX 0

#endif /* LIBSAKURA_LIBSAKURA_CONFIG_H_ */
