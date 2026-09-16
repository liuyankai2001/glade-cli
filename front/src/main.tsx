import { StrictMode } from 'react';
import { createRoot } from 'react-dom/client';
import { Workbench } from './Workbench';
import './style.css';
import './controls.css';

createRoot(document.getElementById('root')!).render(<StrictMode><Workbench /></StrictMode>);
