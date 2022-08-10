/*
 * CustomException.h
 *
 * Copyright © 2019 Deutsches Elektronen-Synchrotron DESY,
 *                       a research centre of the Helmholtz Association.
 *
 * Authors:
 *   2019      Yaroslav Gevorkov <yaroslav.gevorkov@desy.de>
 *
 * This file is part of XGANDALF.
 *
 * XGANDALF is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation, either version 3 of
 * the License, or (at your option) any later version.
 *
 * XGANDALF is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with XGANDALF.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef CUSTOMEXCEPTION_H_
#define CUSTOMEXCEPTION_H_

#include <exception>
#include <string>

namespace xgandalf
{
    //! The Exception class from which every exception type in this project inherits from.
    /*!
     * This is the superclass for all the custom exception types used in this project.
     * This class itself inherits from std::exception.
     */
    class CustomException : public std::exception
    {
      public:
        /*!
         * The constructor taking a message as parameter.
         * @param msg The message of the exception.
         */
        CustomException(const std::string& msg)
            : msg(msg)
        {
        }

        /*!
         * The virtual destructor.
         */
        virtual ~CustomException() throw() {}

        /*!
         * Returns the message of the exception.
         * @return
         */
        virtual const char* what() const throw()
        {
            return msg.c_str();
        }

      private:
        std::string msg;
    };

} // namespace xgandalf
#endif /* CUSTOMEXCEPTION_H_ */
