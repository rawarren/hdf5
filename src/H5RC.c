/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Copyright by The HDF Group.                                               *
 * All rights reserved.                                                      *
 *                                                                           *
 * This file is part of HDF5.  The full HDF5 copyright notice, including     *
 * terms governing use, modification, and redistribution, is contained in    *
 * the files COPYING and Copyright.html.  COPYING can be found at the root   *
 * of the source code distribution tree; Copyright.html can be found at the  *
 * root level of an installed copy of the electronic HDF5 document set and   *
 * is linked from the top-level documents page.  It can also be found at     *
 * http://hdfgroup.org/HDF5/doc/Copyright.html.  If you do not have          *
 * access to either file, you may request a copy from help@hdfgroup.org.     *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/*
 * Programmer:  Mohamad Chaarawi <chaarawi@hdfgroup.org>
 *              August, 2013
 *
 * Purpose:	Read Context APIs to support Exascale FastForward
 *              functionality.
 *              
 */


/****************/
/* Module Setup */
/****************/

/* Interface initialization */
#define H5_INTERFACE_INIT_FUNC	H5RC_init_interface

/***********/
/* Headers */
/***********/
#include "H5private.h"		/* Generic Functions			*/
#include "H5Eprivate.h"		/* Error handling		  	*/
#include "H5ESprivate.h"        /* Event Stacks                         */
#include "H5RCprivate.h"        /* Read Contexts                         */
#include "H5Iprivate.h"		/* IDs			  		*/
#include "H5MMprivate.h"	/* Memory management			*/
#include "H5VLprivate.h"	/* VOL plugins				*/
#include "H5VLiod_client.h"	/* IOD VOL plugin			*/

#ifdef H5_HAVE_EFF

/****************/
/* Local Macros */
/****************/


/******************/
/* Local Typedefs */
/******************/


/********************/
/* Local Prototypes */
/********************/


/*********************/
/* Package Variables */
/*********************/


/*****************************/
/* Library Private Variables */
/*****************************/


/*******************/
/* Local Variables */
/*******************/

/* Declare a free list to manage the H5RC_t struct */
H5FL_DEFINE(H5RC_t);

/* Dataspace ID class */
static const H5I_class_t H5I_RC_CLS[1] = {{
    H5I_RC,        		/* ID class value */
    0,				/* Class flags */
    2,				/* # of reserved IDs for class */
    (H5I_free_t)H5RC_close   	/* Callback routine for closing objects of this class */
}};



/*-------------------------------------------------------------------------
 * Function:	H5RC_init
 *
 * Purpose:	Initialize the interface from some other package.
 *
 * Return:	Success:	non-negative
 *		Failure:	negative
 *
 * Programmer:	Mohamad Chaarawi
 *              August 2013
 *
 *-------------------------------------------------------------------------
 */
herr_t
H5RC_init(void)
{
    herr_t ret_value = SUCCEED;   /* Return value */

    FUNC_ENTER_NOAPI(FAIL)
    /* FUNC_ENTER() does all the work */

done:
    FUNC_LEAVE_NOAPI(ret_value)
} /* end H5RC_init() */


/*--------------------------------------------------------------------------
NAME
   H5RC_init_interface -- Initialize interface-specific information
USAGE
    herr_t H5RC_init_interface()

RETURNS
    Non-negative on success/Negative on failure
DESCRIPTION
    Initializes any interface-specific data or routines.

--------------------------------------------------------------------------*/
static herr_t
H5RC_init_interface(void)
{
    herr_t ret_value = SUCCEED;   /* Return value */

    FUNC_ENTER_NOAPI_NOINIT

    /* Initialize the atom group for the RC IDs */
    if(H5I_register_type(H5I_RC_CLS) < 0)
	HGOTO_ERROR(H5E_DATASPACE, H5E_CANTINIT, FAIL, "unable to initialize interface")

done:
    FUNC_LEAVE_NOAPI(ret_value)
} /* end H5RC_init_interface() */


/*--------------------------------------------------------------------------
 NAME
    H5RC_term_interface
 PURPOSE
    Terminate various H5RC objects
 USAGE
    void H5RC_term_interface()
 RETURNS
    Non-negative on success/Negative on failure
 DESCRIPTION
    Release the atom group and any other resources allocated.
 GLOBAL VARIABLES
 COMMENTS, BUGS, ASSUMPTIONS
     Can't report errors...
 EXAMPLES
 REVISION LOG
--------------------------------------------------------------------------*/
int
H5RC_term_interface(void)
{
    int	n = 0;

    FUNC_ENTER_NOAPI_NOINIT_NOERR

    if(H5_interface_initialize_g) {
	if(H5I_nmembers(H5I_RC) > 0) {
	    (void)H5I_clear_type(H5I_DATASPACE, FALSE, FALSE);
            n++; /*H5I*/
	} /* end if */
        else {
            /* Destroy the dataspace object id group */
	    (void)H5I_dec_type_ref(H5I_RC);
            n++; /*H5I*/

	    /* Shut down interface */
	    H5_interface_initialize_g = 0;
	} /* end else */
    } /* end if */

    FUNC_LEAVE_NOAPI(n)
} /* end H5RC_term_interface() */


/*-------------------------------------------------------------------------
 * Function:	H5Pset_rcapl_version_request
 *
 * Purpose:	Set the request type for acquring the container version, 
 *              whether it has to be exactly the one requested (H5RC_EXACT),
 *              the lowest available greater than the one requested if it is 
 *              not available (H5RC_NEXT), or the last one available ((H5RC_LAST).
 *
 * Return:	Non-negative on success/Negative on failure
 *
 * Programmer:  Mohamad Chaarawi
 *              August 2013
 *
 *-------------------------------------------------------------------------
 */
herr_t
H5Pset_rcapl_version_request(hid_t rcapl_id, H5RC_request_t acquire_req)
{
    H5P_genplist_t *plist;      /* Property list pointer */
    herr_t ret_value = SUCCEED; /* return value */

    FUNC_ENTER_API(FAIL)
    H5TRACE2("e", "iRc", rcapl_id, acquire_req);

    if(H5RC_LAST < acquire_req)
        HGOTO_ERROR (H5E_ARGS, H5E_BADVALUE, FAIL, "Acquire request can be either H5RC_EXACT, H5RC_NEXT, or H5RC_LAST");

    /* Get the plist structure */
    if(NULL == (plist = H5P_object_verify(rcapl_id, H5P_RC_ACQUIRE)))
        HGOTO_ERROR(H5E_ATOM, H5E_BADATOM, FAIL, "can't find object for ID");

    /* Set property */
    if(H5P_set(plist, H5RC_ACQUIRE_CV_REQUEST_NAME, &acquire_req) < 0)
        HGOTO_ERROR(H5E_PLIST, H5E_CANTSET,FAIL, "can't set acquire request");

done:
    FUNC_LEAVE_API(ret_value)
} /* end H5Pset_rcapl_version_request() */


/*-------------------------------------------------------------------------
 * Function:	H5Pget_rcapl_version_request
 *
 * Purpose:	Get the request type for acquring the container version, 
 *
 * Return:	Non-negative on success/Negative on failure
 *
 * Programmer:  Mohamad Chaarawi
 *              September 2013
 *
 *-------------------------------------------------------------------------
 */
herr_t
H5Pget_rcapl_version_request(hid_t rcapl_id, H5RC_request_t *acquire_req)
{
    H5P_genplist_t *plist;      /* Property list pointer */
    herr_t ret_value = SUCCEED; /* return value */

    FUNC_ENTER_API(FAIL)
    H5TRACE2("e", "i*Rc", rcapl_id, acquire_req);

    /* Get the plist structure */
    if(NULL == (plist = H5P_object_verify(rcapl_id, H5P_RC_ACQUIRE)))
        HGOTO_ERROR(H5E_ATOM, H5E_BADATOM, FAIL, "can't find object for ID");

    if(acquire_req) {
        /* Get property */
        if(H5P_get(plist, H5RC_ACQUIRE_CV_REQUEST_NAME, acquire_req) < 0)
            HGOTO_ERROR(H5E_PLIST, H5E_CANTGET,FAIL, "can't get acquire request");
    }

done:
    FUNC_LEAVE_API(ret_value)
} /* end H5Pget_rcapl_version_request() */


/*-------------------------------------------------------------------------
 * Function:	H5RCcreate
 *
 * Purpose:	Wraps an hid_t around a container version and a file ID.
 *              The user is responsible for making sure the container version
 *              is acquired. 
 *
 * Return:	Success:	The ID for a new read context.
 *		Failure:	FAIL
 *
 * Programmer:	Mohamad Chaarawi
 *		August 2013
 *
 *-------------------------------------------------------------------------
 */
hid_t
H5RCcreate(hid_t file_id, uint64_t c_version)
{
    H5VL_object_t *file = NULL;
    H5RC_t *rc = NULL;
    hid_t ret_value;

    FUNC_ENTER_API(FAIL)
    H5TRACE2("i", "iIl", file_id, c_version);

    /* get the file object */
    if(NULL == (file = (H5VL_object_t *)H5I_object_verify(file_id, H5I_FILE)))
        HGOTO_ERROR(H5E_ARGS, H5E_BADVALUE, FAIL, "not a file ID")

    /* create a new read context object */
    if(NULL == (rc = H5RC_create(file, c_version)))
	HGOTO_ERROR(H5E_SYM, H5E_CANTCREATE, FAIL, "unable to create read context")

    rc->vol_cls = file->vol_info->vol_cls;

    /* Get an atom for the RC */
    if((ret_value = H5I_register(H5I_RC, rc, TRUE)) < 0)
	HGOTO_ERROR(H5E_ATOM, H5E_CANTREGISTER, FAIL, "unable to atomize read context handle")

done:
    FUNC_LEAVE_API(ret_value)
} /* end H5RCcreate() */


/*-------------------------------------------------------------------------
 * Function:	H5RC_create
 *
 * Purpose:	Private routine for H5RCcreate
 *
 * Return:	Success:	RC structure
 *		Failure:	NULL
 *
 * Programmer:	Mohamad Chaarawi
 *		August 2013
 *
 *-------------------------------------------------------------------------
 */
H5RC_t *
H5RC_create(void *file, uint64_t c_version)
{
    H5RC_t *rc = NULL;
    H5RC_t *ret_value = NULL;          /* Return value */

    FUNC_ENTER_NOAPI_NOINIT

    /* allocate read context struct */
    if(NULL == (rc = H5FL_CALLOC(H5RC_t)))
	HGOTO_ERROR(H5E_SYM, H5E_NOSPACE, NULL, "can't allocate read context structure")

    rc->file = (H5VL_iod_file_t *)file;
    rc->c_version = c_version;

    rc->req_info.request = NULL;
    rc->req_info.head = NULL;
    rc->req_info.tail = NULL;
    rc->req_info.num_req = 0;

    /* set return value */
    ret_value = rc;

done:
    if(!ret_value && rc)
	rc = H5FL_FREE(H5RC_t, rc);

    FUNC_LEAVE_NOAPI(ret_value)
} /* end H5RC_create() */

herr_t
H5RCget_version(hid_t rc_id, uint64_t *c_version)
{
    H5RC_t *rc = NULL;
    herr_t ret_value = SUCCEED;

    FUNC_ENTER_API(FAIL)
    H5TRACE2("e", "i*Il", rc_id, c_version);

    /* get the RC object */
    if(NULL == (rc = (H5RC_t *)H5I_object_verify(rc_id, H5I_RC)))
        HGOTO_ERROR(H5E_ARGS, H5E_BADVALUE, FAIL, "not a Read Context ID")

    *c_version = rc->c_version;

done:
    FUNC_LEAVE_API(ret_value)
}


/*-------------------------------------------------------------------------
 * Function:	H5RCacquire
 *
 * Purpose:	Acquires a read context from a container for a container
 *              version c_version. The user can specify options to 
 *              what container version to get if the one being asked for
 *              is not available.
 *
 * Return:	Success:	RC ID
 *		Failure:	NULL
 *
 * Programmer:	Mohamad Chaarawi
 *		August 2013
 *
 *-------------------------------------------------------------------------
 */
hid_t
H5RCacquire(hid_t file_id, uint64_t *c_version, 
            hid_t rcapl_id, hid_t estack_id)
{
    H5VL_object_t *file = NULL;
    H5RC_t *rc = NULL;
    H5_priv_request_t  *request = NULL; /* private request struct inserted in event queue */
    void **req = NULL; /* pointer to plugin generate requests (Stays NULL if plugin does not support async */
    hid_t rc_id = FAIL;
    hid_t ret_value = FAIL;

    FUNC_ENTER_API(FAIL)
    H5TRACE4("i", "i*Ilii", file_id, c_version, rcapl_id, estack_id);

    /* get the file object */
    if(NULL == (file = (H5VL_object_t *)H5I_object_verify(file_id, H5I_FILE)))
        HGOTO_ERROR(H5E_ARGS, H5E_BADVALUE, FAIL, "not a file ID")

    /* Get correct property list */
    if(H5P_DEFAULT == rcapl_id)
        rcapl_id = H5P_RC_ACQUIRE_DEFAULT;
    else
        if(TRUE != H5P_isa_class(rcapl_id, H5P_RC_ACQUIRE))
            HGOTO_ERROR(H5E_ARGS, H5E_BADTYPE, FAIL, "not read context acquire property list")

    /* create a new read context object with the provided container
       version. The container version associated with the RC object
       might be changed if the version is not available for read. */
    if(NULL == (rc = H5RC_create(file, *c_version)))
	HGOTO_ERROR(H5E_SYM, H5E_CANTCREATE, FAIL, "unable to create read context")
    rc->vol_cls = file->vol_info->vol_cls;
    /* Get an atom for the RC*/
    if((rc_id = H5I_register(H5I_RC, rc, TRUE)) < 0)
	HGOTO_ERROR(H5E_ATOM, H5E_CANTREGISTER, FAIL, "unable to atomize read context handle")

    if(estack_id != H5_EVENT_STACK_NULL) {
        /* create the private request */
        if(NULL == (request = (H5_priv_request_t *)H5MM_calloc(sizeof(H5_priv_request_t))))
            HGOTO_ERROR(H5E_RESOURCE, H5E_NOSPACE, FAIL, "memory allocation failed");
        req = &request->req;
        request->vol_cls = rc->vol_cls;
    }

    if(H5VL_iod_rc_acquire((H5VL_iod_file_t *)file, rc, c_version, rcapl_id, req) < 0)
        HGOTO_ERROR(H5E_FILE, H5E_CANTINIT, FAIL, "failed to request a read context on container");

    if(request && *req)
        if(H5ES_insert(estack_id, request) < 0)
            HGOTO_ERROR(H5E_SYM, H5E_CANTINIT, FAIL, "failed to insert request in event stack")

    ret_value = rc_id;

done:
    if(ret_value < 0 && rc_id) {
        if(H5I_dec_app_ref(rc_id) < 0)
            HGOTO_ERROR(H5E_SYM, H5E_CANTRELEASE, FAIL, "unable to close read context")
    }
        
    FUNC_LEAVE_API(ret_value)
} /* end H5RCacquire() */


/*-------------------------------------------------------------------------
 * Function:	H5RCrelease
 *
 * Purpose:	Releases a previously acquired container version. Does not
 *              close the RC object (i.e. user has to call H5RCclose()).
 *
 * Return:	Success:	Non-Negative
 *		Failure:	Negative
 *
 * Programmer:	Mohamad Chaarawi
 *		August 2013
 *
 *-------------------------------------------------------------------------
 */
herr_t
H5RCrelease(hid_t rc_id , hid_t estack_id)
{
    H5RC_t *rc = NULL;
    H5_priv_request_t  *request = NULL; /* private request struct inserted in event queue */
    void               **req = NULL; /* pointer to plugin generate requests (Stays NULL if plugin does not support async */
    herr_t ret_value = SUCCEED;

    FUNC_ENTER_API(FAIL)
    H5TRACE2("e", "ii", rc_id, estack_id);

    /* get the RC object */
    if(NULL == (rc = (H5RC_t *)H5I_object_verify(rc_id, H5I_RC)))
        HGOTO_ERROR(H5E_ARGS, H5E_BADVALUE, FAIL, "not a Read Context ID")

    if(estack_id != H5_EVENT_STACK_NULL) {
        /* create the private request */
        if(NULL == (request = (H5_priv_request_t *)H5MM_calloc(sizeof(H5_priv_request_t))))
            HGOTO_ERROR(H5E_RESOURCE, H5E_NOSPACE, FAIL, "memory allocation failed");
        req = &request->req;
        request->vol_cls = rc->vol_cls;
    }

    if(H5VL_iod_rc_release(rc, req) < 0)
        HGOTO_ERROR(H5E_FILE, H5E_CANTINIT, FAIL, "failed to request a release on a read context on container");

    if(request && *req)
        if(H5ES_insert(estack_id, request) < 0)
            HGOTO_ERROR(H5E_SYM, H5E_CANTINIT, FAIL, "failed to insert request in event stack")

done:
    FUNC_LEAVE_API(ret_value)
} /* end H5RCrelease() */


/*-------------------------------------------------------------------------
 * Function:	H5RCpersist
 *
 * Purpose:	Persists a previously acquired container version.
 *
 * Return:	Success:	Non-Negative
 *		Failure:	Negative
 *
 * Programmer:	Mohamad Chaarawi
 *		August 2013
 *
 *-------------------------------------------------------------------------
 */
herr_t
H5RCpersist(hid_t rc_id , hid_t estack_id)
{
    H5RC_t *rc = NULL;
    H5_priv_request_t  *request = NULL; /* private request struct inserted in event queue */
    void               **req = NULL; /* pointer to plugin generate requests (Stays NULL if plugin does not support async */
    herr_t ret_value = SUCCEED;

    FUNC_ENTER_API(FAIL)
    H5TRACE2("e", "ii", rc_id, estack_id);

    /* get the RC object */
    if(NULL == (rc = (H5RC_t *)H5I_object_verify(rc_id, H5I_RC)))
        HGOTO_ERROR(H5E_ARGS, H5E_BADVALUE, FAIL, "not a Read Context ID")

    if(estack_id != H5_EVENT_STACK_NULL) {
        /* create the private request */
        if(NULL == (request = (H5_priv_request_t *)H5MM_calloc(sizeof(H5_priv_request_t))))
            HGOTO_ERROR(H5E_RESOURCE, H5E_NOSPACE, FAIL, "memory allocation failed");
        req = &request->req;
        request->vol_cls = rc->vol_cls;
    }

    if(H5VL_iod_rc_persist(rc, req) < 0)
        HGOTO_ERROR(H5E_FILE, H5E_CANTINIT, FAIL, "failed to request a persist on a read context on container");

    if(request && *req)
        if(H5ES_insert(estack_id, request) < 0)
            HGOTO_ERROR(H5E_SYM, H5E_CANTINIT, FAIL, "failed to insert request in event stack")

done:
    FUNC_LEAVE_API(ret_value)
} /* end H5RCpersist() */


/*-------------------------------------------------------------------------
 * Function:	H5RCsnapshot
 *
 * Purpose:	Creates a permanently visible snapshot of a container state.
 *
 * Return:	Success:	Non-Negative
 *		Failure:	Negative
 *
 * Programmer:	Mohamad Chaarawi
 *		August 2013
 *
 *-------------------------------------------------------------------------
 */
herr_t
H5RCsnapshot(hid_t rc_id , const char *snapshot_name, hid_t estack_id)
{
    H5RC_t *rc = NULL;
    H5_priv_request_t  *request = NULL; /* private request struct inserted in event queue */
    void               **req = NULL; /* pointer to plugin generate requests (Stays NULL if plugin does not support async */
    herr_t ret_value = SUCCEED;

    FUNC_ENTER_API(FAIL)
    H5TRACE3("e", "i*si", rc_id, snapshot_name, estack_id);

    /* Check arguments */
    if(!snapshot_name || !*snapshot_name)
	HGOTO_ERROR(H5E_ARGS, H5E_BADVALUE, FAIL, "invalid snapshot name")

    /* get the RC object */
    if(NULL == (rc = (H5RC_t *)H5I_object_verify(rc_id, H5I_RC)))
        HGOTO_ERROR(H5E_ARGS, H5E_BADVALUE, FAIL, "not a Read Context ID")

    if(estack_id != H5_EVENT_STACK_NULL) {
        /* create the private request */
        if(NULL == (request = (H5_priv_request_t *)H5MM_calloc(sizeof(H5_priv_request_t))))
            HGOTO_ERROR(H5E_RESOURCE, H5E_NOSPACE, FAIL, "memory allocation failed");
        req = &request->req;
        request->vol_cls = rc->vol_cls;
    }

    if(H5VL_iod_rc_snapshot(rc, snapshot_name, req) < 0)
        HGOTO_ERROR(H5E_FILE, H5E_CANTINIT, FAIL, "failed to request a snapshot on a read context on container");

    if(request && *req)
        if(H5ES_insert(estack_id, request) < 0)
            HGOTO_ERROR(H5E_SYM, H5E_CANTINIT, FAIL, "failed to insert request in event stack")

done:
    FUNC_LEAVE_API(ret_value)
} /* end H5RCsnapshot() */


/*-------------------------------------------------------------------------
 * Function:	H5RCclose
 *
 * Purpose:	Closes the specified read context ID.  The ID will no longer be
 *		valid for accessing the read context.
 *
 * Return:	Non-negative on success/Negative on failure
 *
 * Programmer:	Mohamad Chaarawi
 *		August 2013
 *
 *-------------------------------------------------------------------------
 */
herr_t
H5RCclose(hid_t rc_id)
{
    herr_t ret_value = SUCCEED;         /* Return value */

    FUNC_ENTER_API(FAIL)
    H5TRACE1("e", "i", rc_id);

    /* Check args */
    if(NULL == H5I_object_verify(rc_id,H5I_RC))
	HGOTO_ERROR(H5E_ARGS, H5E_BADTYPE, FAIL, "not an read context ID")

    if(H5I_dec_app_ref(rc_id) < 0)
    	HGOTO_ERROR(H5E_SYM, H5E_CANTRELEASE, FAIL, "unable to close read context")

done:
    FUNC_LEAVE_API(ret_value)
} /* end H5RCclose() */


/*-------------------------------------------------------------------------
 * Function:	H5RC_close
 *
 * Purpose:	Closes the specified read context ID.  The ID will no longer be
 *		valid for accessing the read context.
 *
 * Return:	Non-negative on success/Negative on failure
 *
 * Programmer:	Mohamad Chaarawi
 *		August 2013
 *
 *-------------------------------------------------------------------------
 */
herr_t
H5RC_close(H5RC_t *rc)
{
    FUNC_ENTER_NOAPI_NOINIT_NOERR

    if(rc->req_info.num_req) {
        H5VL_iod_request_t *cur_req = rc->req_info.head;
        H5VL_iod_request_t *next_req = NULL;
        H5VL_iod_request_t *prev;
        H5VL_iod_request_t *next;

        while(cur_req) {
            next_req = cur_req->trans_next;

            /* remove the current request from the linked list */
            prev = cur_req->trans_prev;
            next = cur_req->trans_next;
            if (prev) {
                if (next) {
                    prev->trans_next = next;
                    next->trans_prev = prev;
                }
                else {
                    prev->trans_next = NULL;
                    rc->req_info.tail = prev;
                }
            }
            else {
                if (next) {
                    next->trans_prev = NULL;
                    rc->req_info.head = next;
                }
                else {
                    rc->req_info.head = NULL;
                    rc->req_info.tail = NULL;
                }
            }

            cur_req->trans_prev = NULL;
            cur_req->trans_next = NULL;

            rc->req_info.num_req --;

            H5VL_iod_request_decr_rc(cur_req);

            cur_req = next_req;
        }
        HDassert(0 == rc->req_info.num_req);
    }

    rc = H5FL_FREE(H5RC_t, rc);

    FUNC_LEAVE_NOAPI(SUCCEED)
} /* end H5RC_close() */

#endif /* H5_HAVE_EFF */
