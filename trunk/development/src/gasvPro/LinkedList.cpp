/**
 * Copyright 2012 Benjamin Raphael, Suzanne Sindi, Anthony Cannistra, Hsin-Ta Wu, Luke Peng, Selim Onal
 *
 *  This file is part of the GASVPro code distribution.
 * 
 *  GASVPro is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 * 
 *  GASVPro is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *  
 *  You should have received a copy of the GNU General Public License
 *  along with gasv.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

#include "LinkedList.h"

        /**
         * Construct the list.
         */
        template <class Object>
        List<Object>::List( )
        {
            header = new ListNode<Object>;
       }

        /**
         * Copy constructor.
         */
        template <class Object>
        List<Object>::List( const List<Object> & rhs )
        {
            header = new ListNode<Object>;
            *this = rhs;
        }

        /**
         * Destructor.
         */
        template <class Object>
        List<Object>::~List( )
        {
            makeEmpty( );
            delete header;
        }

        /**
         * Test if the list is logically empty.
         * Return true if empty, false, otherwise.
         */
        template <class Object>
        bool List<Object>::isEmpty( ) const
        {
            return header->next == NULL;
        }

        /**
         * Make the list logically empty.
         */
        template <class Object>
        void List<Object>::makeEmpty( )
        {
            while( !isEmpty( ) )
                remove( first( ).retrieve( ) );
        }

        /**
         * Return an iterator representing the header node.
         */
        template <class Object>
        ListItr<Object> List<Object>::zeroth( ) const
        {
            return ListItr<Object>( header );
        }

        /**
         * Return an iterator representing the first node in the list.
         * This operation is valid for empty lists.
         */
        template <class Object>
        ListItr<Object> List<Object>::first( ) const
        {
            return ListItr<Object>( header->next );
        }

        /**
         * Insert item x after p.
         */
        template <class Object>
        void List<Object>::insert( const Object & x, const ListItr<Object> & p )
        {
            if( p.current != NULL )
                p.current->next = new ListNode<Object>( x, p.current->next );
        }

        /**
         * Return iterator corresponding to the first node containing an item x.
         * Iterator isPastEnd if item is not found.
         */
        template <class Object>
        ListItr<Object> List<Object>::find( const Object & x ) const
        {
/* 1*/      ListNode<Object> *itr = header->next;

/* 2*/      while( itr != NULL && itr->element != x )
/* 3*/          itr = itr->next;

/* 4*/      return ListItr<Object>( itr );
        }

        /**
         * Return iterator prior to the first node containing an item x.
         */
        template <class Object>
        ListItr<Object> List<Object>::findPrevious( const Object & x ) const
        {
/* 1*/      ListNode<Object> *itr = header;

/* 2*/      while( itr->next != NULL && itr->next->element != x )
/* 3*/          itr = itr->next;

/* 4*/      return ListItr<Object>( itr );
        }

        /**
         * Remove the first occurrence of an item x.
         */
        template <class Object>
        void List<Object>::remove( const Object & x )
        {
            ListItr<Object> p = findPrevious( x );

            if( p.current->next != NULL )
            {
                ListNode<Object> *oldNode = p.current->next;
                p.current->next = p.current->next->next;  // Bypass deleted node
                delete oldNode;
            }
        }

        /**
         * Deep copy of linked lists.
         */
        template <class Object>
        const List<Object> & List<Object>::operator=( const List<Object> & rhs )
        {
            if( this != &rhs )
            {
                makeEmpty( );

                ListItr<Object> ritr = rhs.first( );
                ListItr<Object> itr = zeroth( );
                for( ; !ritr.isPastEnd( ); ritr.advance( ), itr.advance( ) )
                    insert( ritr.retrieve( ), itr );
            }
            return *this;
        }

