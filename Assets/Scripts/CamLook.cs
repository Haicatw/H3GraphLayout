using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class CamLook : MonoBehaviour
{
    float horizontalSpeed = 2.0f;
    float verticalSpeed = 2.0f;
    float speed = 3.0f;

    Camera cam;

    void Start()
    {
        cam = GetComponent<Camera>();
        //freeze Z
    }

    // Update is called once per frame
    void Update()
    {
        mouseMovement();
        mouseSelect();
    }

    void mouseMovement()
    {
        // Get the mouse delta. This is not in the range -1...1
        float h = horizontalSpeed * Input.GetAxis("Mouse X");
        float v = verticalSpeed * -Input.GetAxis("Mouse Y");

        transform.Rotate(v, h, 0);
    }

    void mouseSelect()
    {
        if (Input.GetMouseButtonDown(0))
        {
            RaycastHit hit;
            Ray ray = cam.ScreenPointToRay(Input.mousePosition);
            if (Physics.Raycast(ray, out hit, 100.0f))
            {
                Debug.Log("Select: " + hit.transform.gameObject.name);
                GameObject graph = GameObject.FindGameObjectWithTag("GlobalManager");
                GraphLayout graphScript = graph.GetComponent<GraphLayout>();
                graphScript.updateTranslation(hit.transform.gameObject.name);
            }
        }
    }
}