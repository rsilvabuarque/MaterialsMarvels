'use client';

import { useEffect, useState } from 'react'
import styles from "./page.module.css";

export default function Home() {
  const [message, setMessage] = useState("");
  const [loading, setLoading] = useState(true);

  fetch('/api/')
      .then(res => res.json())
      .then(data => {
          setMessage(data.hello);
          setLoading(false);
      })

    return (
        <div className={styles.container}>
            <p> {!loading ? message : "Loading.."}</p>
        </div>
    );
}
